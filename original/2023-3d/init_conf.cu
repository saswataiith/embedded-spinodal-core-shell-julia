#include"binary.h"

__global__ void  kernel_real2Complex(cufftDoubleComplex *d_comp, cufftDoubleComplex *d_phi,
				     cufftDoubleComplex *d_dfdc, cufftDoubleComplex *d_dfdphi,
 				     double *d_realc, double *d_realphi, int nx, int ny, int nz){
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
   unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
   unsigned int tid = iz + nz * (iy + ix * ny);

   if (tid < nx * ny * nz) {
     d_comp[tid].x = d_realc[tid];
     d_phi[tid].x = d_realphi[tid];
     d_comp[tid].y = 0.0;
     d_phi[tid].y = 0.0;
   
     d_dfdc[tid].x = d_comp[tid].x;
     d_dfdc[tid].y = d_comp[tid].y;
     d_dfdphi[tid].x = d_phi[tid].x;
     d_dfdphi[tid].y = d_phi[tid].y;
   }
}
__global__ void desiredRandom_no(double *d_random, 
				double noise, int nx, int ny, int nz){  
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
   unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
   unsigned int tid = iz + nz * (iy + ix * ny);
  
   if (tid < nx * ny * nz) 
      d_random[tid] = (2.0 * d_random[tid] - 1.0) * noise;
     
}

void Init_conf(){

  curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen, SEED);

  curandGenerateUniformDouble(gen, d_randm, nx * ny * nz);
  
  dim3 threads (threads_x, threads_y, threads_z);
 
  dim3 blocks(blocks_x, blocks_y, blocks_z);

  printf("Kernel configuration<<<(%d, %d, %d), (%d, %d, %d)>>>\n", 
		  blocks.x, blocks.y, blocks.z, threads.x, threads.y, threads.z);

  CHECK (cudaDeviceSynchronize());
  desiredRandom_no<<<blocks, threads>>>(d_randm, noise, nx, ny, nz);
  CHECK (cudaDeviceSynchronize());
  CHECK (cudaGetLastError());

  CHECK (cudaMemcpy(randm, d_randm, nBytes, cudaMemcpyDeviceToHost));

  double sumrnd, rr;

  sumrnd = 0.0;
  counter = 0;
  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      for(int k = 0; k < nz; k++){
       rr = ((double) i - (double) (nx_half )) * ((double) i - (double) (nx_half)) 
           + ((double) j - (double) (ny_half )) * ((double) j - (double) (ny_half))
           + ((double) k - (double) (nz_half )) * ((double) k - (double) (nz_half));
       rr = sqrt(rr);

       if (rr < R)
        phi[k + nz * (j + i * ny)].x = 1.0;

        phi[k + nz * (j + i * ny)].y = 0.0;

       if (phi[k + nz * (j + i * ny)].x >= phicutoff){
         sumrnd += randm[k + nz * (j + i * ny)];
         counter++; 
       }
      }
    }
  }

  mean = sumrnd / (double) counter;
  printf("mean of random number inside the particle = %le\n", mean);

  one_by_counter = 1.0 / (double)counter;

  CHECK (cudaMemcpy(d_dfdphi, phi, nBytesComplex, cudaMemcpyHostToDevice));
  CHECK (cudaMemcpy(d_phi, phi, nBytesComplex, cudaMemcpyHostToDevice));

  sum = 0.0;
  sumrnd = 0.0;

  int countbeta1 = 0;
  int countbeta2 = 0;

  for(int i = 0; i < nx; i++){
   for(int j = 0; j < ny; j++){
    for(int k = 0; k < nz; k++){
     if (phi[k + nz * (j + i * ny)].x > phicutoff){
       
       if (flag == 1)
        comp[k + nz * (j + i * ny)].x = c_0 * (c_beta1 + c_beta2) + randm[k + nz * (j + i * ny)] - mean;

       if (flag == 2){
  	  if (k <= nz_half){
            comp[k + nz * (j + i * ny)].x = c_beta1;
	    countbeta1++;
	  }else{
            comp[k + nz * (j + i * ny)].x = c_beta2;
	    countbeta2++;
	  }
       }

	 sumrnd += (randm[k + nz * (j + i * ny)] - mean);
         sum += comp[k + nz * (j + i * ny)].x;
     } else
         comp[k + nz * (j + i * ny)].x = c_alpha;

       comp[k + nz * (j + i * ny)].y = 0.0;
    }
   }
  }

  mean = sum * one_by_counter; 
  printf("Particle avg comp: %le\n", mean);
  mean = sumrnd * one_by_counter;
  printf("mean of (random number - mean) inside the particle: %le\n", mean);
  if (flag == 2) printf("Initally beta1 points: %d and beta2 points: %d\n", countbeta1, countbeta2);

  CHECK (cudaMemcpy(d_dfdc, comp, nBytesComplex, cudaMemcpyHostToDevice));
  CHECK (cudaMemcpy(d_comp, comp, nBytesComplex, cudaMemcpyHostToDevice));
   
  // Show result 
  FILE *fpq;
  fpq = fopen("XY_profile.in","w");
  int k = nz_half;
  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
     fprintf(fpq, "%d\t%d\t%le\t%le\n", i, j, comp[k + nz * (j + i * ny)].x, phi[k + nz * (j + i * ny)].x);  
    } fprintf(fpq,"\n");
  }
  fclose(fpq);
  fpq = fopen("ZY_profile.in","w");
  int i = nx_half;
    for(int k = 0; k < nz; k++){
  for(int j = 0; j < ny; j++){
     fprintf(fpq, "%d\t%d\t%le\t%le\n", j, k, comp[k + nz * (j + i * ny)].x, phi[k + nz * (j + i * ny)].x);  
    } fprintf(fpq,"\n");
  }
  fclose(fpq);
  fpq = fopen("ZX_profile.in","w");
  int j = ny_half;
  for(int i = 0; i < nx; i++){
    for(int k = 0; k < nz; k++){
     fprintf(fpq, "%d\t%d\t%le\t%le\n", i, k, comp[k + nz * (j + i * ny)].x, phi[k + nz * (j + i * ny)].x);  
    } fprintf(fpq,"\n");
  }
  fclose(fpq);

   CHECK (cudaFree(d_randm));
   free(randm);
}

void Read_Restart(){
  FILE *fpread;
  char fr[100];
  
  double *realphi, *realc, *d_realphi, *d_realc;
  realphi = (double*)malloc(nBytes);
  realc = (double*)malloc(nBytes);
  CHECK (cudaMalloc((void **)&d_realc, nBytes));
  CHECK (cudaMalloc((void **)&d_realphi, nBytes));
  
  sprintf (fr,"conf.%09d", initcount);
  fpread = fopen (fr, "rb");
  if(fread (&realc[0], sizeof(double), nx * ny * nz, fpread));
  if(fread (&realphi[0], sizeof(double), nx * ny * nz, fpread));
  fclose (fpread);
 
  /*
  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      for(int k = 0; k < nz; k++){
       comp[k + nz * (j + i * ny)].x = realc[k + nz * (j + i * ny)];
       phi[k + nz * (j + i * ny)].x = realphi[k + nz * (j + i * ny)];
       comp[k + nz * (j + i * ny)].y = 0.0;
       phi[k + nz * (j + i * ny)].y = 0.0;
      }
    }
  }
  */

  CHECK (cudaMemcpy(d_realc, realc, nBytes, cudaMemcpyHostToDevice));
  CHECK (cudaMemcpy(d_realphi, realphi, nBytes, cudaMemcpyHostToDevice));

  dim3 threads (threads_x, threads_y, threads_z);
 
  dim3 blocks(blocks_x, blocks_y, blocks_z);

  kernel_real2Complex<<<blocks, threads>>>(d_comp, d_phi, d_dfdc, d_dfdphi, d_realc, d_realphi, nx, ny, nz);

  free(realc);
  free(realphi);
  cudaFree(d_realc);
  cudaFree(d_realphi);
}
