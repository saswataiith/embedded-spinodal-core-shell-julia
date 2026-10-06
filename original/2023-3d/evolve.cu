
#include"binary.h"
#include"EvoKernels.h"

void Evolve(){
  void Output_conf(int count);

  double   *tempreal, *d_tempreal, *tempphi, *d_tempphi;
  double   alloycomp, total, err; 
  int	   *boundsFlag, *d_boundsFlag;
  cufftDoubleComplex *d_mu;
  void     *t_storage = NULL;
  size_t   t_storage_bytes = 0;
  double   maxerror, *maxerr_d;
  double   *sum_d;
  
  CHECK(cudaMalloc(&t_storage, t_storage_bytes));
  CHECK(cudaMalloc((void**)&maxerr_d, sizeof(double)));
  CHECK(cudaMalloc((void**)&sum_d, sizeof(double)));
  CHECK (cudaMalloc((void **)&d_mu, nBytesComplex));
  CHECK (cudaMalloc((void **)&d_tempreal, nBytes));
  
  cub::DeviceReduce::Max(t_storage, t_storage_bytes, d_tempreal, maxerr_d, nx * ny * nz);
  
  cudaMalloc(&t_storage, t_storage_bytes);

  boundsFlag = (int *)malloc(sizeof(int));

  *boundsFlag = 1;

  CHECK (cudaMalloc((void **)&d_boundsFlag, sizeof(int)));
  tempreal = (double *) malloc(nBytes);
  //tempphi = (double *) malloc(nBytes);

  sum = 0.0;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
       tempreal[k + nz * (j + i * ny)] = comp[k + nz * (j + i * ny)].x;
    //   tempphi[k + nz * (j + i * ny)] = phi[k + nz * (j + i * ny)].x;
       sum = sum + comp[k + nz * (j + i * ny)].x;
      }
    }
  }

  alloycomp = sum * one_by_nxnynz;

  CHECK (cudaMemcpy(d_tempreal, tempreal,  nBytes, cudaMemcpyHostToDevice));
  //CHECK (cudaMemcpy(d_tempphi, tempphi,  nBytes, cudaMemcpyHostToDevice));
  free(phi);
  free(tempreal);
  //free(comp);

  cufftHandle plan = 0;
  CHECK_CUFFT (cufftPlan3d(&plan, nx, ny, nz, CUFFT_Z2Z));

  CHECK (cudaMemcpy(d_boundsFlag, boundsFlag, sizeof(int), cudaMemcpyHostToDevice));

  dim3 threads (threads_x, threads_y, threads_z);
 
  dim3 blocks(blocks_x, blocks_y, blocks_z);

  printf("Kernel configuration<<<(%d, %d, %d), (%d, %d, %d)>>>\n", 
		  blocks.x, blocks.y, blocks.z, threads.x, threads.y, threads.z);

  int loop_condition = 1;

  CHECK_CUFFT (cufftExecZ2Z(plan, d_phi, d_phi, CUFFT_FORWARD));
    CHECK_CUFFT (cufftExecZ2Z(plan, d_comp, d_comp, CUFFT_FORWARD));

  double iStart, tElaps;

  iStart = cpuSecond();

  for(int count = 0; count <= total_steps; count++){

    if( ((count % print_steps_1 == 0) && (count <= total_steps_1) ) ||
        ((count % print_steps_2 == 0) && (count <= total_steps_2) ) ||
        ((count % print_steps) == 0) || (count == total_steps) 
	||   loop_condition == 0 || *boundsFlag == 0 ){
       printf("Total_time = %le\n", sim_time);
       printf("Writing files!. Time = %d\n", count);
       Output_conf(count);
    }

    if (loop_condition == 0 || *boundsFlag == 0 || count > total_steps)
       break;

    kernel_cphi_der<<<blocks, threads>>>(d_dfdc, d_dfdphi, mobility,
          	                         A, A0, B, chi, P, c_alpha, 
        				 c_beta1, c_beta2, nx, ny, nz);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    //Execute a complex-to-complex 3d FFT in-plane transformation
    CHECK_CUFFT (cufftExecZ2Z(plan, d_dfdc, d_dfdc, CUFFT_FORWARD));

    Compute_mu<<<blocks, threads>>>(d_mu, d_comp, d_dfdc, d_kx2, d_ky2, d_kz2, kappa_c, nx, ny, nz);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    CHECK_CUFFT (cufftExecZ2Z(plan, d_mu, d_mu, CUFFT_INVERSE));

    Muscaling<<<blocks, threads>>>(d_mu, one_by_nxnynz, nx, ny, nz);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    ComputeGradmu_x<<<blocks, threads>>>(d_gradmu, d_mu, one_by_nxnynz, nx, ny, nz, dx);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

     ComputeGrad_x<<<blocks, threads>>>(d_dfdc, d_gradmu, d_dfdphi, nx, ny, nz, dx);
     CHECK (cudaDeviceSynchronize());
     CHECK (cudaGetLastError());

    ComputeGradmu_y<<<blocks, threads>>>(d_gradmu, d_mu, one_by_nxnynz, nx, ny, nz, dy);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

     ComputeGrad_y<<<blocks, threads>>>(d_dfdc, d_gradmu, d_dfdphi, nx, ny, nz, dy);
     CHECK (cudaDeviceSynchronize());
     CHECK (cudaGetLastError());

    ComputeGradmu_z<<<blocks, threads>>>(d_gradmu, d_mu, one_by_nxnynz, nx, ny, nz, dy);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

     ComputeGrad_z<<<blocks, threads>>>(d_dfdc, d_gradmu, d_dfdphi, nx, ny, nz, dz);
     CHECK (cudaDeviceSynchronize());
     CHECK (cudaGetLastError());

    reset_dfdphi<<<blocks, threads>>>(d_dfdphi, nx, ny, nz);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    CHECK_CUFFT (cufftExecZ2Z(plan, d_dfdphi, d_dfdphi, CUFFT_FORWARD));
    CHECK_CUFFT (cufftExecZ2Z(plan, d_dfdc, d_dfdc, CUFFT_FORWARD));

    kernel_comp_phi<<<blocks, threads>>>(d_comp, d_dfdc, 
       	                                 d_phi, d_dfdphi, d_kx1, d_ky1, d_kz1, 
					 d_kx2, d_ky2, d_kz2, dt, kappa_c, 
       			     		 kappa_phi, Q, relax_coeff, nx, ny, nz);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());
  
    CHECK (cudaMemcpy(comp, d_dfdc, nBytesComplex, cudaMemcpyDeviceToHost));
    total = comp[0].x * one_by_nxnynz;
    err = fabs(total - alloycomp);
    if(err > COMPERR){
       printf("Elements NOT conserved...... Exiting\n");
       printf("error = %le\n", err);
       exit(0);
    }

    CHECK_CUFFT (cufftExecZ2Z(plan, d_dfdc, d_dfdc, CUFFT_INVERSE));
    CHECK_CUFFT (cufftExecZ2Z(plan, d_dfdphi, d_dfdphi, CUFFT_INVERSE));

    Normalization<<<blocks, threads>>>(d_dfdc, d_dfdphi, one_by_nxnynz, nx, ny, nz);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    Bounds<<<blocks, threads>>>(d_dfdc, d_dfdphi, d_boundsFlag, nx, ny, nz);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    CHECK (cudaMemcpy(boundsFlag, d_boundsFlag, sizeof(int), cudaMemcpyDeviceToHost));
    if (*boundsFlag == 0) printf("Composition is out of bounds. Exited\n");
   
    MaxError<<<blocks, threads>>>(d_tempreal, d_dfdc, nx, ny, nz);
    cub::DeviceReduce::Max(t_storage, t_storage_bytes, d_tempreal, maxerr_d, nx * ny * nz);       
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    CHECK(cudaMemcpy(&maxerror, maxerr_d, sizeof(double), cudaMemcpyDeviceToHost));
    //printf("Maxerror = %le\n", maxerror); 
    
    if(maxerror <= Tolerance){
      printf("!!!! CONVERGENCE ACHIEVED !!!!.\n"); 
      printf("Maxerror at convergence = %le\n", maxerror); 
      loop_condition = 0;      
    }
     
    UpdateTempreal<<<blocks, threads>>>(d_tempreal, d_dfdc, nx, ny, nz);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

   sim_time = sim_time + dt;
   
  }

  cudaFree(d_mu);
  free(comp);
  //cudaFree(d_tempphi);
  //cudaFree(d_tempreal);
  tElaps = cpuSecond() - iStart;
  printf("Time consumed in 'count' loop = %lf seconds\n", tElaps);
}
