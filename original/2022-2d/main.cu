
#include"binary.h"
#include"get_input.cu"
#include"init_conf.cu"
#include"evolve.cu"
#include"out_conf.cu"

int main(){

   void Get_Input_Parameters(char *fnin, char *fnout);	
   void Init_conf();
   void Read_Restart();
   void kvector(double *kx1, double *ky1, double *kx2, double *ky2);
   void Evolve();
   void devPropQuery();

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

   devPropQuery();

   nx_half = nx/2;
   ny_half = ny/2;

   // nBytesComplex is a size_t type of cuDoubleComplex variable which defines the size of the mem required
   nBytesComplex = nx * ny * sizeof(cufftDoubleComplex);
   // nBytes is a size_t type of double variable which defines the size of the mem required
   nBytes = sizeof(double) * nx * ny;
   
   one_by_nxny = 1.0 /((double)nx * ny);

   // Allocate memory on host //
   randm = (double *)malloc(nBytes);
   comp = (cufftDoubleComplex *)malloc(nBytesComplex);
   
   // Allocate memory on device //
   CHECK (cudaMalloc((void **)&d_randm, nBytes));
   CHECK (cudaMalloc((void **)&d_comp, nBytesComplex));

   //Below 3 vairables work in both CPU and GPU memory. 
   //Not required to copy their values into CPU variable to print.
   cudaMallocManaged(&phi, nBytesComplex);
   cudaMallocManaged(&d_dfdphi, nBytesComplex);
   cudaMallocManaged(&d_dfdc, nBytesComplex);
  
   CHECK_CUFFT (cufftPlan2d(&plan, nx, ny, CUFFT_Z2Z));

   if (initflag == 0){
      Init_conf();  
   }else{
      Read_Restart();
   }

  
   kx1 = (double *) malloc(sizeof(double) * nx);
   ky1 = (double *) malloc(sizeof(double) * ny);
   kx2 = (double *) malloc(sizeof(double) * nx);
   ky2 = (double *) malloc(sizeof(double) * ny);

   CHECK (cudaMalloc((void **)&d_kx1, sizeof(double) * nx));
   CHECK (cudaMalloc((void **)&d_ky1, sizeof(double) * ny));
   CHECK (cudaMalloc((void **)&d_kx2, sizeof(double) * nx));
   CHECK (cudaMalloc((void **)&d_ky2, sizeof(double) * ny));

   kvector(kx1, ky1, kx2, ky2);

   CHECK (cudaMemcpy(d_kx1, kx1, sizeof(double) * nx, cudaMemcpyHostToDevice));
   CHECK (cudaMemcpy(d_ky1, ky1, sizeof(double) * ny, cudaMemcpyHostToDevice));
   CHECK (cudaMemcpy(d_kx2, kx2, sizeof(double) * nx, cudaMemcpyHostToDevice));
   CHECK (cudaMemcpy(d_ky2, ky2, sizeof(double) * ny, cudaMemcpyHostToDevice));

   Evolve();
   
   fclose(fpout);

   /* Cleanup */
   if (initflag == 0)
   CURAND_CALL (curandDestroyGenerator(gen));
   CHECK(cudaFree(d_comp));
   CHECK (cudaFree(d_dfdc));
   CHECK (cudaFree(d_dfdphi));
   CHECK (cudaFree(d_randm));
   free(comp);
   free(randm);
   free(kx1);
   free(ky1);
   free(kx2);
   free(ky2);
   CHECK (cudaFree(d_kx1));
   CHECK (cudaFree(d_ky1));
   CHECK (cudaFree(d_kx2));
   CHECK (cudaFree(d_ky2));

   //Reset device
   cudaDeviceReset();
   return EXIT_SUCCESS;

}

void kvector(double *kx1, double *ky1, double *kx2, double *ky2){
  double dkx, dky; 

  dkx = 2.0 * M_PI/((double)nx * dx);
  dky = 2.0 * M_PI/((double)ny * dy);

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

}

void devPropQuery(){
  
   cudaDeviceProp deviceProp;
   CHECK (cudaGetDeviceProperties(&deviceProp, dev));
   printf(" Using device %d: %s,\n", dev, deviceProp.name);
   CHECK (cudaSetDevice(dev));
   if (!deviceProp.concurrentManagedAccess) UMflag = 0;

   if (!UMflag) printf("Unified memory is not supported in the chosen GPU\n");

}
