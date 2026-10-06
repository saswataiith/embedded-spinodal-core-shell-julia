#include"binary.h"

__global__ void RealToComplex(cufftDoubleComplex *d_comp, double *c,
	                      cufftDoubleComplex *d_dfdc,	
		              cufftDoubleComplex *d_dfdphi, double *phi1,
			      int nx, int ny){
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
   unsigned int tid = iy + ix * ny;
  
   if (tid <= nx * ny) {
    d_comp[tid].x = c[tid];
    d_dfdc[tid].x = c[tid];
    d_dfdphi[tid].x = phi1[tid];
    d_comp[tid].y = 0.0; 
    d_dfdc[tid].y = 0.0;
    d_dfdphi[tid].y = 0.0;
   }
}

__global__ void desiredRandom_no(double *d_random, 
				double noise, int nx, int ny){  
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
   unsigned int tid = iy + ix * ny;
  
   if (tid <= nx * ny) 
      d_random[tid] = (2.0 * d_random[tid] - 1.0) * noise;
     
}

void Init_conf(){

  curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen, SEED);

  curandGenerateUniformDouble(gen, d_randm, nx * ny);
  
  dim3 threads (threads_x, threads_y);  // no of threads per block
  dim3 blocks( (nx+threads.x-1)/threads.x, 
		(ny+threads.y-1)/threads.y); // no of blocks per grid
  printf("Kernel configuration<<<(%d, %d), (%d, %d)>>>\n", 
		blocks.x, blocks.y, threads.x, threads.y);

  CHECK (cudaDeviceSynchronize());
  desiredRandom_no<<<blocks, threads>>>(d_randm, noise, nx, ny);
  CHECK (cudaDeviceSynchronize());
  CHECK (cudaGetLastError());

  CHECK (cudaMemcpy(randm, d_randm, nBytes, cudaMemcpyDeviceToHost));

  double sumrnd, rr, rr1, rr2;
  double a = 0.5, b = 1.0;
  double j1, j2;
  double ch_height, ch_width;
  double width, width1;

  if (flag == 1){
   ch_width = (double)nx;
   ch_height = (R * R)/ch_width;
   j1 = ((double)ny - ch_height) / 2.0;
   j2 = j1 + ch_height;

   j1 = ceil(j1);
   j2 = ceil(j2);
   printf("j1 %lf, j2 %lf\n", j1, j2);
   width1 = 2.0 * pre_factor;
  }
  
  width = 2.0 * pre_factor;
  printf("width = %lf\n", width);

  sumrnd = 0.0;
  counter = 0;
  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
     if (flag == 0){ 
      rr = ((double) i - (double) (nx_half )) 
          * ((double) i - (double) (nx_half))
         + ((double) j - (double) (ny_half ))
         * ((double) j - (double) (ny_half ));
      rr = sqrt(rr);
      phi[j + i * ny].x = a * (b - tanh((rr - R) / width));

      // phi[j + i * ny].x = 1.0;
     }

     if(flag == 1){
      phi[j + i * ny].x = a * ( b + tanh(( (double)j - j1) / width1) ) - 
			  a * ( b + tanh(( (double)j - j2) / width1) );
     }
    
     if(flag == 2){
      rr =  ( ((double) i - (double) (nx_half)) / (major/minor) )
          * ( ((double) i - (double) (nx_half)) / (major/minor) )
          + ( (double) j - (double) (ny_half))
          * ( (double) j - (double) (ny_half));
      rr = sqrt(rr);
      phi[j + i * ny].x = a * (b - tanh( ( rr - (R/sqrt(major/minor)) ) / width));
     }

     if(flag == 3){
      rr =   ((double) i - (double) (nx_half)) 
           * ((double) i - (double) (nx_half))
           + ((double) j - (double) (ny_half))
           * ((double) j - (double) (ny_half));
      rr = sqrt(rr);
      phi[j + i * ny].x = a * (b - tanh((rr - Ro) / width)) - 
		          a * (b - tanh((rr - R)/width));
     }

     if(flag == 4){
       rr = (double)(i - nx_half) * (double)(i - nx_half) * (double)(i - nx_half) * (double)(i - nx_half)  +
            (double)(j - ny_half) * (double)(j - ny_half) * (double)(j - ny_half) * (double)(j - ny_half) ;
       rr = pow(rr, 0.25);

       phi[j + i * ny].x = 0.5 * (1.0 - tanh((rr - R) / width));
     }

     if(flag == 5){
      rr1 =  ( ((double) i - (double) (nx_half)) / (major/minor) ) * dx
          * ( ((double) i - (double) (nx_half)) / (major/minor) ) * dx
          + ( (double) j - (double) (ny_half)) * dy
          * ( (double) j - (double) (ny_half)) * dy;
      rr1 = sqrt(rr1);
      rr2 = ((double) i - (double) (nx_half )) * dx
          * ((double) i - (double) (nx_half)) * dx
         + ((double) j - (double) (ny_half )) * dy
         * ((double) j - (double) (ny_half )) * dy;
      rr2 = sqrt(rr2);
      phi[j + i * ny].x = a * (b - tanh( ( rr1 - (Ro/sqrt(major/minor)) ) / width)) - a * (b - tanh((rr2 - R)/width));
     }

     if(flag == 6){
      rr1 =  ( ((double) i - (double) (nx_half)) / (major/minor) ) 
          * ( ((double) i - (double) (nx_half)) / (major/minor) ) 
          + ( (double) j - (double) (ny_half)) 
          * ( (double) j - (double) (ny_half)) ;
      rr1 = sqrt(rr1);
      rr2 = ((double) i - (double) (nx_half ))
          * ((double) i - (double) (nx_half)) 
         + ((double) j - (double) (ny_half ))
         * ((double) j - (double) (ny_half ));
      rr2 = sqrt(rr2);
      phi[j + i * ny].x = - a * (b - tanh( ( rr1 - (R/sqrt(major/minor)) ) / width)) + 
	      		    a * (b - tanh((rr2 - Ro)/width));
     }

      if (phi[j + i * ny].x >= phicutoff){
         sumrnd += randm[j + i * ny];
         counter++; 
      }
        
         phi[j + i * ny].y = 0.0;
    }
  }
  
  mean = sumrnd / (double) counter;
  printf("mean of random number inside the particle = %le\n", mean);

  one_by_counter = 1.0 / (double)counter;

  printf("counter = %d\n", counter);

  CHECK (cudaMemcpy(d_dfdphi, phi, nBytesComplex, cudaMemcpyHostToDevice));

  sum = 0.0;
  sumrnd = 0.0;
 
  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      if (phi[j + i * ny].x >= phicutoff){
        comp[j + i * ny].x = c_0 * (c_beta1 + c_beta2) + randm[j + i * ny] - mean;
	sumrnd += (randm[j + i * ny] - mean);
        sum += comp[j + i * ny].x;
      } else 
        comp[j + i * ny].x = c_alpha;

      comp[j + i * ny].y = 0.0;
    }
  }
 
  mean = sum * one_by_counter; 
  printf("Particle avg comp: %le\n", mean);
  mean = sumrnd * one_by_counter;
  printf("mean of (random number - mean) inside the particle: %le\n", mean);

  CHECK (cudaMemcpy(d_dfdc, comp, nBytesComplex, cudaMemcpyHostToDevice));
  CHECK (cudaMemcpy(d_comp, comp, nBytesComplex, cudaMemcpyHostToDevice));
   
  // Show result 
  FILE *fpq;
  fpq = fopen("profile.in","w");
  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
     fprintf(fpq, "%le\t%le\t%le\t%le\n", i * dx, j * dy, comp[j + i * ny].x, phi[j + i * ny].x);  
    } fprintf(fpq,"\n");
  }
  fclose(fpq);
  
 char fn[100];
  sprintf(fn, "profile.lcut_%0.3lf",dx);
  fpq = fopen(fn, "w");
  for(int j = 0; j < ny; j++)
    fprintf(fpq,"%le\t%le\t%lf\n", j * dy, comp[j + nx/2 * ny].x, phi[j + nx/2 * ny].x);
  fclose(fpq);

  CHECK (cudaFree(phi));

}

void Read_Restart()
{
  FILE *fpread;
  char fr[100];
/*
  int dimx = 1;
  int dimy = 512;
  
  dim3 threads (dimx, dimy);
  dim3 blocks( (nx + threads.x-1) / threads.x, (ny + threads.y-1)/threads.y);
 
  double *d_c, *d_phi1, *c, *phi1;
  CHECK (cudaMalloc((void **)&d_c, nBytes));
  CHECK (cudaMalloc((void **)&d_phi1, nBytes));
  c = (double*)malloc(sizeof(double) * nx * ny); 
  phi1 = (double*)malloc(sizeof(double) * nx * ny); 

  sprintf (fr,"conf.%09d", initcount);
  fpread = fopen (fr, "r");
  if(fread (&c[0], sizeof(double), nx * ny, fpread));
  if(fread (&phi1[0], sizeof(double), nx * ny, fpread));
  fclose (fpread);
  
  CHECK (cudaMemcpy(d_phi1, phi1, nBytes, cudaMemcpyHostToDevice));
  CHECK (cudaMemcpy(d_c, c, nBytes, cudaMemcpyHostToDevice));

  RealToComplex<<<blocks, threads>>>(d_comp, d_c, d_dfdc, d_dfdphi, d_phi1);

  counter = 0;
  for(int i = 0; i < nx; i++)
    for(int j = 0; j < ny; j++)
      if (phi1[j + i * ny] >= phicutoff)
            counter++;

  one_by_counter = 1.0 / (double)counter;

  CHECK (cudaFree(phi));
  CHECK (cudaFree(d_phi1));
  CHECK (cudaFree(d_c));
  free(c);
  free(phi1);
  */

  sprintf (fr,"conf.%11d", initcount);
  fpread = fopen (fr, "r");
  if(fread (&comp[0], sizeof(double), 2 * nx * ny, fpread));
  if(fread (&phi[0], sizeof(double), 2 *nx * ny, fpread));
  fclose (fpread);
  
  counter = 0;
  for(int i = 0; i < nx; i++)
    for(int j = 0; j < ny; j++)
      if (phi[j + i * ny].x >= phicutoff)
            counter++;

  one_by_counter = 1.0 / (double)counter;
  
  CHECK (cudaMemcpy(d_comp, comp, nBytesComplex, cudaMemcpyHostToDevice));
  CHECK (cudaMemcpy(d_dfdc, comp, nBytesComplex, cudaMemcpyHostToDevice));
  CHECK (cudaMemcpy(d_dfdphi, phi, nBytesComplex, cudaMemcpyHostToDevice));
  CHECK (cudaFree(phi));
}
