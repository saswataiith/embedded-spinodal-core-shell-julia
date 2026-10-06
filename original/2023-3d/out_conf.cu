
#include"binary.h"

__global__ void  kernel_Complex2real(cufftDoubleComplex *d_dfdc, 
				     cufftDoubleComplex *d_dfdphi,
 				     double *d_realc, double *d_mux, 
					int nx, int ny, int nz){

   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
   unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
   unsigned int tid = iz + nz * (iy + ix * ny);

   if (tid < nx * ny * nz) {
     d_realc[tid] = d_dfdc[tid].x;
     d_mux[tid] = d_dfdphi[tid].x;
   }
}

void Output_conf(int steps){
 
  FILE *q;
  char fn[100];

  double *realphi, *realc, *d_realc;
  realphi = (double*)malloc(nBytes);
  realc = (double*)malloc(nBytes);
  CHECK (cudaMalloc((void **)&d_realc, nBytes));
  
  dim3 threads (threads_x, threads_y, threads_z);
 
  dim3 blocks(blocks_x, blocks_y, blocks_z);

  kernel_Complex2real<<<blocks, threads>>>(d_dfdc, d_dfdphi, d_realc, d_gradmu, nx, ny, nz);

  CHECK (cudaMemcpy(realc, d_realc, nBytes, cudaMemcpyDeviceToHost));
  CHECK (cudaMemcpy(realphi, d_gradmu, nBytes, cudaMemcpyDeviceToHost));

  sprintf(fn, "conf.%09d", steps + initcount);
  q = fopen(fn, "wb");
  fwrite (&realc[0], sizeof(double), nx * ny * nz, q);
  fwrite (&realphi[0], sizeof(double), nx * ny * nz , q);
  fclose(q);

  free(realc);
  free(realphi);
  cudaFree(d_realc);
}
