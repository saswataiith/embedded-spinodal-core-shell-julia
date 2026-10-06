
#include"binary.h"

__global__ void ComplexToReal(cufftDoubleComplex *d_comp,
		              double *d_c,
		              cufftDoubleComplex *d_dfdphi,
			      double *d_phi1, int nx, int ny){

   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
   unsigned int tid = iy + ix * ny;

   if (tid < nx * ny) {
        d_c[tid] = d_comp[tid].x;
	d_phi1[tid] = d_dfdphi[tid].x;
   }
}

void Output_conf(int steps){
 
  FILE *q;
  char fn[100];
  //int ind;
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

  ComplexToReal<<<blocks, threads>>>(d_comp, d_c, d_dfdphi, d_phi1);

  CHECK (cudaMemcpy(phi1, d_phi1, nBytes, cudaMemcpyDeviceToHost));
  CHECK (cudaMemcpy(c, d_c, nBytes, cudaMemcpyDeviceToHost));
*/
  CHECK (cudaMemcpy(comp, d_comp, nBytesComplex, cudaMemcpyDeviceToHost));
  //if (steps % 100000 == 0){
  //sprintf(fn, "/home/pankaj/Src_1/R40/conf.%07d", steps + initcount);
  sprintf(fn, "conf.%09d", steps + initcount);
  q = fopen(fn, "w");
  fwrite (&comp[0], sizeof(double), 2 * nx * ny, q);
  fwrite (&d_dfdphi[0], sizeof(double), 2 * nx * ny , q);
  fclose(q);
  //}

  //sprintf(fn, "prof_gp.%09d", steps + initcount);
  //q = fopen(fn, "w");
  //for (int i = 0; i < nx; i++){
  // for (int j = 0; j < ny; j++){
  //    ind = j + i * ny;
  //    fprintf(q, "%d\t%d\t%le\t%le\n", i, j, d_dfdc[ind].x, d_dfdphi[ind].x);		 
  // } fprintf(q, "\n");
  //} fclose(q);

  //CHECK (cudaFree(d_phi1));
  //CHECK (cudaFree(d_c));
  //free(c);
  //free(phi1);
}
