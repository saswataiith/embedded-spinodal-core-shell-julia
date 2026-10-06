
#include"binary.h"

__global__ void Random_no(double *d_random, 
			  cufftDoubleComplex *d_comp,
				int ny, int nx){  
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
   unsigned int tid = iy + ix * ny;

   if (tid < nx * ny) {
//      d_random[tid] =  d_comp[tid].x * (1.0 - d_comp[tid].x) 
//             * (2.0 * d_random[tid] - 1.0);
      d_random[tid] = (2.0 * d_random[tid] - 1.0);
   }
     
}

__global__ void dfdc_kernel(cufftDoubleComplex *d_dfdc, 
		cufftDoubleComplex *d_comp,
                double *d_hphi, double *d_gphi,	
		double A, double B, double chi, double P, 
		double c_alpha, double c_beta1, double c_beta2,
		int nx, int ny){     
     double  ctemp;
     unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
     unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
     unsigned int tid = iy + ix * ny;
  
  if(tid < nx * ny){
     ctemp = d_comp[tid].x;

     d_dfdc[tid].x = 2.0 * A * (1.0 - d_hphi[tid]) * (ctemp - c_alpha) + 
	             2.0 * B * d_hphi[tid] * (ctemp - c_beta1) * (ctemp - c_beta2) *
		     (2.0 * ctemp - c_beta1 - c_beta2) - chi * P * d_gphi[tid];
     d_dfdc[tid].y = 0.0;

  }
}
    
__global__ void derivative_mu(cufftDoubleComplex *d_comp, cufftDoubleComplex *d_dfdc,  
		cufftDoubleComplex * d_mux, cufftDoubleComplex * d_muy, 
       	        double *d_kx1, double *d_ky1, double *d_kx2, double *d_ky2, 
		double kappa_c, int nx, int ny){

  double kpow2, rc, fp, mu; 
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int tid = iy + ix * ny;
 
  if(tid < nx * ny){
    kpow2 = d_kx2[ix] * d_kx2[ix] + d_ky2[iy] * d_ky2[iy];
    rc = d_comp[tid].x;
    fp = d_dfdc[tid].x;
    mu = fp + 2.0 * kappa_c * kpow2 * rc;
    d_mux[tid].y = 1.0 * d_kx1[ix] * mu;
    d_muy[tid].y = 1.0 * d_ky1[iy] * mu;

    rc = d_comp[tid].y;
    fp = d_dfdc[tid].y;
    mu = fp + 2.0 * kappa_c * kpow2 * rc;
    d_mux[tid].x = - 1.0 * d_kx1[ix] * mu;
    d_muy[tid].x = - 1.0 * d_ky1[iy] * mu;
  }
}

__global__ void mu_scaling(cufftDoubleComplex *d_mux, cufftDoubleComplex *d_muy, 
		cufftDoubleComplex *d_dfdphi, double *d_hphi, double mobility, double mob_tol,
	  	double  one_by_nxny, int nx, int ny){
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int tid = iy + ix * ny;
  
  //double temp_phi, hphi1;

  if(tid < nx * ny){
    d_mux[tid].x *= one_by_nxny;
    d_muy[tid].x *= one_by_nxny;
    d_mux[tid].y = 0.0;
    d_muy[tid].y = 0.0;
/*
    if (d_dfdphi[tid].x < 0.05){
     d_mux[tid].x *= (mobility * d_hphi[tid] + mob_tol);
     d_muy[tid].x *= (mobility * d_hphi[tid] + mob_tol);
    }else{
     d_mux[tid].x *= mobility * d_hphi[tid];
     d_muy[tid].x *= mobility * d_hphi[tid];
    }
    */
     d_mux[tid].x *= mobility * d_hphi[tid];
     d_muy[tid].x *= mobility * d_hphi[tid];

     //temp_phi = d_dfdphi[tid].y;
     //hphi1 = temp_phi * temp_phi * temp_phi * (10.0 - 15.0 * temp_phi + 6.0 * temp_phi * temp_phi);
     d_mux[tid].y = 0.0;//*= mobility * hphi1;
     d_muy[tid].y = 0.0;//*= mobility * hphi1;
/*
     temp_phi = d_dfdphi[tid].y;
     hphi1 = temp_phi * temp_phi * temp_phi * (10.0 - 15.0 * temp_phi + 6.0 * temp_phi * temp_phi);
    if (d_dfdphi[tid].y < 0.05){
     d_mux[tid].y *= 0.0;//mobility * hphi1;// + mob_tol);
     d_muy[tid].y *= 0.0;//mobility * hphi1;// + mob_tol);
    } else {
     d_mux[tid].y *= 0.0;//mobility * hphi1;// + mob_tol);
     d_muy[tid].y *= 0.0;//mobility * hphi1;// + mob_tol);

    }
    */
  }

}

__global__ void CH_kernel(cufftDoubleComplex *d_comp, cufftDoubleComplex *d_dfdc, 
		cufftDoubleComplex *d_mux, cufftDoubleComplex *d_muy,
       	        double *d_kx1, double *d_ky1, double *d_kx2, double *d_ky2, 
		double dt, double kappa_c, double Q, int nx, int ny){
  double abcxIm, abcxRe, abcyIm, abcyRe,abc; 
  double rc, fp, rhs, lhs, kpow2, kpow4;
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int tid = iy + ix * ny;

  if(tid < nx * ny){
    kpow2 = d_kx2[ix] * d_kx2[ix] + d_ky2[iy] * d_ky2[iy];
    kpow4 = kpow2 * kpow2;
    lhs = 1.0 + 2.0 * Q * kappa_c * kpow4 * dt;

    abcxIm = 1.0 * d_kx1[ix] * d_mux[tid].x;
    abcxRe = - 1.0 * d_kx1[ix] * d_mux[tid].y;

    abcyIm = 1.0 * d_ky1[iy] * d_muy[tid].x;
    abcyRe = - 1.0 * d_ky1[iy] * d_muy[tid].y;

    abc = abcxRe + abcyRe;
    rc = d_comp[tid].x;
    fp = dt * abc + (1.0 + 2.0 * Q * kappa_c * kpow4 * dt) * rc;
    rhs = fp;
    d_comp[tid].x = rhs / lhs;

    abc = abcxIm + abcyIm;
    rc = d_comp[tid].y;
    fp = dt * abc + (1.0 + 2.0 * Q * kappa_c * kpow4 * dt) * rc;
    rhs = fp;
    d_comp[tid].y = rhs / lhs;
  }
}

__global__ void Normalization(cufftDoubleComplex *d_x,
                         double one_by_nxny, int nx, int ny){
     unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
     unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
     unsigned int tid = iy + ix * ny;
     
  if(tid < nx * ny){
     d_x[tid].x *= one_by_nxny;
     d_x[tid].y = 0.0;
  }
}

__global__ void Bounds(cufftDoubleComplex *d_comp, 
			int *d_boundsFlag, int nx, int ny){
     unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
     unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
     unsigned int tid = iy + ix * ny;

  if(tid < nx * ny){
    if(d_comp[tid].x < -0.4 || d_comp[tid].x > 1.4)
      *d_boundsFlag = 0;
  }
}

__global__ void MaxError(double *dtempreal, cufftDoubleComplex *d_comp,
			int nx, int ny){
     unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
     unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
     unsigned int tid = iy + ix * ny;

  if(tid < nx * ny)
    dtempreal[tid] = fabs(dtempreal[tid] - d_comp[tid].x);

}

__global__ void UpdateTempreal(double *dtempreal, cufftDoubleComplex *d_comp,
			       int nx, int ny){
     unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
     unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
     unsigned int tid = iy + ix * ny;

  if(tid < nx * ny)
    dtempreal[tid] = d_comp[tid].x;

}
      
__global__ void AddNoise(cufftDoubleComplex *d_comp, cufftDoubleComplex *d_dfdphi, 
				double *d_randm, double mean, 
				double noise, double phicutoff, int nx, int ny){
     unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
     unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
     unsigned int tid = iy + ix * ny;
  double temp;

  if(tid < nx * ny){
    if (d_dfdphi[tid].x >= phicutoff){
      temp = d_randm[tid] - mean;
      d_comp[tid].x = d_comp[tid].x + noise * temp;
    } else
      d_comp[tid].x = d_comp[tid].x ;
  }
}

__global__ void gradphi_x_y(cufftDoubleComplex *d_dfdphi, 
			cufftDoubleComplex *d_gradphix,
			cufftDoubleComplex *d_gradphiy,
			double *d_kx1, double *d_ky1, int nx, int ny){
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
   unsigned int tid = iy + ix * ny;
 
   if (tid < nx * ny) {
     d_gradphix[tid].x = -1.0 * d_kx1[ix] * d_dfdphi[tid].y;
     d_gradphix[tid].y =  d_kx1[ix] * d_dfdphi[tid].x;


     d_gradphiy[tid].x = -1.0 * d_ky1[iy] * d_dfdphi[tid].y;
     d_gradphiy[tid].y =  d_ky1[iy] * d_dfdphi[tid].x;

   }	   
}

__global__ void phiGrad(cufftDoubleComplex *d_gradphi, 
			cufftDoubleComplex *d_gradphix, 
			cufftDoubleComplex *d_gradphiy, 
			int nx, int ny){
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
   unsigned int tid = iy + ix * ny;

   if (tid < nx * ny) {
     d_gradphi[tid].x = sqrt(d_gradphix[tid].x * d_gradphix[tid].x 
				+ d_gradphiy[tid].x * d_gradphiy[tid].x);
     d_gradphi[tid].y = sqrt(d_gradphix[tid].y * d_gradphix[tid].y 
				+ d_gradphiy[tid].y * d_gradphiy[tid].y);

   }	   
}

void Evolve(){
  void Output_conf(int count);
  void Convergence(cufftDoubleComplex *d_dfdc, double *tempreal, 
		  int count, int *loop_condition);
  
  double   *tempreal, *dtempreal;
  double   alloycomp, total, err; 
  int	   *boundsFlag, *d_boundsFlag;
  cufftDoubleComplex *d_mux, *d_muy;
  double   *hphi, *gphi, *d_hphi, *d_gphi, ptemp;
  void     *t_storage = NULL;
  size_t   t_storage_bytes = 0;
  double   maxerror, *maxerr_d;
  double   *sum_d;
  double   *meanVar1 = (double *)malloc(nBytes); 
  cufftDoubleComplex *d_gradphix, *d_gradphiy, *d_gradphi, *gradphi;
  
  gradphi = (cufftDoubleComplex *)malloc(nBytesComplex);
  CHECK(cudaMalloc(&t_storage, t_storage_bytes));
  CHECK(cudaMalloc((void**)&maxerr_d, sizeof(double)));
  CHECK(cudaMalloc((void**)&sum_d, sizeof(double)));
  CHECK (cudaMalloc((void **)&d_gradphix, nBytesComplex));
  CHECK (cudaMalloc((void **)&d_gradphiy, nBytesComplex));
  CHECK (cudaMalloc((void **)&d_gradphi, nBytesComplex));
  CHECK (cudaMalloc((void **)&d_mux, nBytesComplex));
  CHECK (cudaMalloc((void **)&d_muy, nBytesComplex));
  CHECK (cudaMalloc((void **)&d_hphi, nBytes));
  CHECK (cudaMalloc((void **)&d_gphi, nBytes));
  CHECK (cudaMalloc((void **)&dtempreal, nBytes));

  cub::DeviceReduce::Max(t_storage, t_storage_bytes, dtempreal, maxerr_d, nx * ny);
  
  cudaMalloc(&t_storage, t_storage_bytes);

  boundsFlag = (int *)malloc(sizeof(int));

  *boundsFlag = 1;

  if (initflag == 1){
   curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
   curandSetPseudoRandomGeneratorSeed(gen, SEED);
  }

  tempreal = (double *) malloc(nBytes);
  hphi = (double *) malloc(nBytes);
  gphi = (double *) malloc(nBytes);

  CHECK (cudaMalloc((void **)&d_boundsFlag, sizeof(int)));

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      tempreal[j + i * ny] = comp[j + i * ny].x;
      ptemp = d_dfdphi[j + i * ny].x;
      hphi[j + i * ny] = ptemp * ptemp * ptemp * (10.0 - 15.0 * ptemp + 6.0 * ptemp * ptemp);
      gphi[j + i * ny] = (ptemp * ptemp ) * (1.0 - ptemp) * (1.0 - ptemp);
    }
  }
  sum = 0.0;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      sum = sum + d_dfdc[j + i * ny].x;
    }
  }

  alloycomp = sum * one_by_nxny;

  printf("Particle avg comp after summing Fourier modes = %le\n", alloycomp);

  CHECK (cudaMemcpy(d_hphi, hphi,  nBytes, cudaMemcpyHostToDevice));
  CHECK (cudaMemcpy(d_gphi, gphi,  nBytes, cudaMemcpyHostToDevice));
  CHECK (cudaMemcpy(dtempreal, tempreal,  nBytes, cudaMemcpyHostToDevice));

  free(hphi);
  free(gphi);
  free(tempreal);

  CHECK (cudaMemcpy(d_boundsFlag, boundsFlag, sizeof(int), cudaMemcpyHostToDevice));

  dim3 threads (threads_x, threads_y);
  dim3 blocks( (nx + threads.x-1) / threads.x, (ny + threads.y-1)/threads.y);

  int loop_condition = 1;

  double iStart, tElaps;

  iStart = cpuSecond();
//  double inputnoise = noise;

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

    dfdc_kernel<<<blocks, threads>>>(d_dfdc, d_comp, d_hphi, d_gphi,
       	                         A, B, chi, P, c_alpha, 
       				 c_beta1, c_beta2, nx, ny);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    //Execute a complex-to-complex 2d FFT in-plane transformation
    CHECK_CUFFT (cufftExecZ2Z(plan, d_dfdc, d_dfdc, CUFFT_FORWARD));
    CHECK_CUFFT (cufftExecZ2Z(plan, d_comp, d_comp, CUFFT_FORWARD));

    derivative_mu<<<blocks, threads>>>(d_comp, d_dfdc, d_mux, d_muy, 
        	     d_kx1, d_ky1, d_kx2, d_ky2, kappa_c, nx, ny);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    CHECK_CUFFT (cufftExecZ2Z(plan, d_mux, d_mux, CUFFT_INVERSE));
    CHECK_CUFFT (cufftExecZ2Z(plan, d_muy, d_muy, CUFFT_INVERSE));

    mu_scaling<<<blocks, threads>>>(d_mux, d_muy, d_dfdphi, d_hphi, mobility, 
		    	             mob_tol, one_by_nxny, nx, ny);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    CHECK_CUFFT (cufftExecZ2Z(plan, d_mux, d_mux, CUFFT_FORWARD));
    CHECK_CUFFT (cufftExecZ2Z(plan, d_muy, d_muy, CUFFT_FORWARD));

    CH_kernel<<<blocks, threads>>>(d_comp, d_dfdc, d_mux, d_muy,
       	                     d_kx1, d_ky1, d_kx2, d_ky2, dt, kappa_c, 
       			     Q, nx, ny);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    //Conservation of mass
    CHECK (cudaMemcpy(comp, d_comp, nBytesComplex, cudaMemcpyDeviceToHost));
    total = comp[0].x * one_by_nxny;
    err = fabs(total - alloycomp);
    if(err > COMPERR){
       printf("Elements NOT conserved...... Exiting\n");
       printf("error = %le\n", err);
       exit(0);
    }

    CHECK_CUFFT (cufftExecZ2Z(plan, d_comp, d_comp, CUFFT_INVERSE));

    Normalization<<<blocks, threads>>>(d_comp, one_by_nxny, nx, ny);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

     Bounds<<<blocks, threads>>>(d_comp, d_boundsFlag, nx, ny);
     CHECK (cudaDeviceSynchronize());
     CHECK (cudaGetLastError());

     CHECK (cudaMemcpy(boundsFlag, d_boundsFlag, sizeof(int), cudaMemcpyDeviceToHost));
     if (*boundsFlag == 0) printf("Composition is out of bounds. Exited\n");

    //Convergence check in GPU using cub library
    MaxError<<<blocks, threads>>>(dtempreal, d_comp, nx, ny);
    cub::DeviceReduce::Max(t_storage, t_storage_bytes, dtempreal,
                         maxerr_d, nx*ny);       
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

    CHECK(cudaMemcpy(&maxerror, maxerr_d, sizeof(double), cudaMemcpyDeviceToHost));
    //printf("Maxerror = %le\n", maxerror); 
    
    if(maxerror <= Tolerance){
      printf("!!!! CONVERGENCE ACHIEVED !!!!.\n"); 
      printf("Maxerror at convergence = %le\n", maxerror); 
      loop_condition = 0;      
    }

/*
    if ( (count < noise_steps && initflag == 0) || (count > noise_steps && (count % NatEveryStp) == 0)){
//      if (count < 1000) noise = 0.01;
      //else noise = inputnoise;

      curandGenerateUniformDouble(gen, d_randm, nx * ny);

      Random_no<<<blocks, threads>>>(d_randm, d_comp, nx, ny);
      CHECK (cudaDeviceSynchronize());
      CHECK (cudaGetLastError());
     
      CHECK(cudaMemcpy(meanVar1, d_randm, nx * ny * sizeof(double), cudaMemcpyDeviceToHost));

      sum = 0.0;
      for(int i=0; i<nx; i++)
        for(int j=0; j<ny; j++)
          if (d_dfdphi[j + i * ny].x >= phicutoff)
            sum += meanVar1[j + i * ny];

      mean =  sum * one_by_counter; 

      AddNoise<<<blocks, threads>>>(d_comp, d_dfdphi, d_randm, mean, noise, 
		               phicutoff , nx, ny);
      CHECK (cudaDeviceSynchronize());
      CHECK (cudaGetLastError());
*/
/*
      sum = 0.0;
      for(int i=0; i<nx; i++)
        for(int j=0; j<ny; j++)
          if (d_dfdphi[j + i * ny].x >= phicutoff)
            sum += (meanVar1[j + i * ny] - mean);

      mean =  sum * one_by_counter; 
      printf("mean of (randm - mean) %le\n", mean);

      CHECK (cudaMemcpy(comp, d_comp, nBytesComplex, cudaMemcpyDeviceToHost));

      FILE *fp;
      fp = fopen("Noisecheck","w");
      for (int i=0; i<nx; i++){
         for (int j=0; j<ny; j++){
           fprintf(fp,"%d\t%d\t%le\n", i, j, comp[j+i*ny].x);
	 } fprintf(fp,"\n");
      }fclose(fp);
      */
//  }

   UpdateTempreal<<<blocks, threads>>>(dtempreal, d_comp, nx, ny);
    CHECK (cudaDeviceSynchronize());
    CHECK (cudaGetLastError());

   sim_time = sim_time + dt;
   
  }

  cudaFree(d_hphi);
  cudaFree(d_gphi);
  cudaFree(d_mux);
  cudaFree(d_gradphix);
  cudaFree(d_gradphiy);
  cudaFree(d_muy);
  cudaFree(dtempreal);
  cudaFree(d_gradphi);
  free(gradphi);
  tElaps = cpuSecond() - iStart;
  printf("Time consumed in 'count' loop = %lf seconds\n", tElaps);
}
