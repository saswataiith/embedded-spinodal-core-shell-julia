#include"binary.h"
#include"get_input.cu"
#include"init_conf.cu"
#include"evolve.cu"
#include"out_conf.cu"

int main(){

   void Get_Input_Parameters(char *fnin, char *fnout);	
   void Init_conf();
   void Read_Restart();
   void kvector(double *kx1, double *ky1, double *kz1, 
		   double *kx2, double *ky2, double *kz2);
   void Evolve();

//   system("rm -rf /home/pankaj/Src_1/R40/conf*");

   char finput[15] = "bin1ary";
   char fnin[30], fnout[30];

   FILE *fp;

   if (!(fp = fopen(finput, "r"))) {
     printf("File: %s could not opened\n", finput);
     exit(EXIT_FAILURE);
   }

   if (fscanf(fp, "%s", fnin) == 1) {
     printf("Input Parameters Filename: %s\n", fnin);
   }

   if (fscanf(fp, "%s", fnout) == 1) {
     printf("Output Parameters Filename: %s\n", fnout);
   }

   if (!(fpout = fopen (fnout, "w"))) {
     printf ("File:%s could not be opened\n", fnout);
     exit (EXIT_FAILURE);
   }

   fclose(fp);

   Get_Input_Parameters(fnin, fnout);
   
   //setup GPU device.
   //int dev = 1;
   //cudaDeviceProp deviceProp;
   //CHECK (cudaGetDeviceProperties(&deviceProp, dev));
   //printf(" Using device %d: %s,\n", dev, deviceProp.name);
   //CHECK (cudaSetDevice(dev));

   nx_half = nx/2;
   ny_half = ny/2;
   nz_half = nz/2;

   // nBytesComplex is a size_t type of cuDoubleComplex variable which defines the size of the mem required
   nBytesComplex = nx * ny * nz * sizeof(cufftDoubleComplex);
   // nBytes is a size_t type of double variable which defines the size of the mem required
   nBytes = sizeof(double) * nx * ny * nz;

   one_by_nxnynz = 1.0 /((double)nx * ny * nz);

   // Allocate memory on host //
   randm = (double *)malloc(nBytes);
   comp = (cufftDoubleComplex *)malloc(nBytesComplex);
   phi = (cufftDoubleComplex *)malloc(nBytesComplex);
   
   // Allocate memory on device //
   CHECK (cudaMalloc((void **)&d_randm, nBytes));
   CHECK (cudaMalloc((void **)&d_comp, nBytesComplex));
   CHECK (cudaMalloc((void **)&d_dfdc, nBytesComplex));
   CHECK (cudaMalloc((void **)&d_phi, nBytesComplex));
   CHECK (cudaMalloc((void **)&d_dfdphi, nBytesComplex));
   CHECK (cudaMalloc((void **)&d_gradmu, nBytes));

   blocks_x = ceil(nx/threads_x);
   blocks_y = ceil(ny/threads_y);
   blocks_z = ceil(nz/threads_z);

   if (initflag == 0){
      Init_conf();  
   }else{
      Read_Restart();
   }
  
   kx1 = (double *) malloc(sizeof(double) * nx);
   ky1 = (double *) malloc(sizeof(double) * ny);
   kz1 = (double *) malloc(sizeof(double) * nz);
   kx2 = (double *) malloc(sizeof(double) * nx);
   ky2 = (double *) malloc(sizeof(double) * ny);
   kz2 = (double *) malloc(sizeof(double) * nz);

   CHECK (cudaMalloc((void **)&d_kx1, sizeof(double) * nx));
   CHECK (cudaMalloc((void **)&d_ky1, sizeof(double) * ny));
   CHECK (cudaMalloc((void **)&d_kz1, sizeof(double) * nz));
   CHECK (cudaMalloc((void **)&d_kx2, sizeof(double) * nx));
   CHECK (cudaMalloc((void **)&d_ky2, sizeof(double) * ny));
   CHECK (cudaMalloc((void **)&d_kz2, sizeof(double) * nz));

   kvector(kx1, ky1, kz1, kx2, ky2, kz2);

   CHECK (cudaMemcpy(d_kx1, kx1, sizeof(double) * nx, cudaMemcpyHostToDevice));
   CHECK (cudaMemcpy(d_ky1, ky1, sizeof(double) * ny, cudaMemcpyHostToDevice));
   CHECK (cudaMemcpy(d_kz1, kz1, sizeof(double) * nz, cudaMemcpyHostToDevice));
   CHECK (cudaMemcpy(d_kx2, kx2, sizeof(double) * nx, cudaMemcpyHostToDevice));
   CHECK (cudaMemcpy(d_ky2, ky2, sizeof(double) * ny, cudaMemcpyHostToDevice));
   CHECK (cudaMemcpy(d_kz2, kz2, sizeof(double) * nz, cudaMemcpyHostToDevice));

   Evolve();
   
   fclose(fpout);

   // Cleanup 
   CHECK(cudaFree(d_comp));
   CHECK (cudaFree(d_dfdc));
   CHECK (cudaFree(d_dfdphi));
   CHECK (cudaFree(d_phi));
   CHECK (cudaFree(d_gradmu));
   //free(comp); //freed in evolve.cu if not checking for conservation
   free(kx1);
   free(ky1);
   free(kz1);
   free(kx2);
   free(ky2);
   free(kz2);
   CHECK (cudaFree(d_kx1));
   CHECK (cudaFree(d_ky1));
   CHECK (cudaFree(d_kz1));
   CHECK (cudaFree(d_kx2));
   CHECK (cudaFree(d_ky2));
   CHECK (cudaFree(d_kz2));

   //Reset device
   cudaDeviceReset();
   return EXIT_SUCCESS;

}

void kvector(double *kx1, double *ky1, double *kz1,
		double *kx2, double *ky2, double *kz2){
  double dkx, dky, dkz; 

  dkx = 2.0 * M_PI/((double)nx * dx);
  dky = 2.0 * M_PI/((double)ny * dy);
  dkz = 2.0 * M_PI/((double)nz * dz);

  for (int i = 0; i < nx; i++) {
   if (i < nx_half){
      kx1[i] = (double) i * dkx;
      kx2[i] = (double) i * dkx;
   }
   else if (i == nx_half) {
      kx1[i] = 0.0;
      kx2[i] = (double) i * dkx;
   }
   else {
      kx1[i] = (double) (i - nx) * dkx;
      kx2[i] = (double) (i - nx) * dkx;
   }
  }

  for (int j = 0; j < ny; j++) {
   if (j < ny_half){
      ky1[j] = (double) j * dky;
      ky2[j] = (double) j * dky;
   }
   else if (j == ny_half) {
      ky1[j] = 0.0;
      ky2[j] = (double) j * dky;
   }
   else {
      ky1[j] = (double) (j - ny) * dky;
      ky2[j] = (double) (j - ny) * dky;
   }
  }

  for (int k = 0; k < nz; k++) {
   if (k < nz_half){
      kz1[k] = (double) k * dkz;
      kz2[k] = (double) k * dkz;
   }
   else if (k == nz_half) {
      kz1[k] = 0.0;
      kz2[k] = (double) k * dkz;
   }
   else {
      kz1[k] = (double) (k - nz) * dkz;
      kz2[k] = (double) (k - nz) * dkz;
   }
  }
}
