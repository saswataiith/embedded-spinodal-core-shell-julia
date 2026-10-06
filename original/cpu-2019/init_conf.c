#include"binary.h"

void Init_Conf()
{
  FILE *fp;
  char fn[100];

  double *random_num;
  double sum, mean;

  random_num = (double*) malloc(sizeof(double) * nx * ny);

  /* Phi profile */
  for (int i = 0; i < nx; i++){ 
    for (int j = 0; j < ny; j++) {
       if ( j>= nx/3 && j<2*nx/3)
//       if((double)(i - nx_half)* (double)(i - nx_half) + (double)(j - ny_half) * (double)(j - ny_half) <= R * R )
         phi[j + i * ny][Re] = 1.0;
       else 
         phi[j + i * ny][Re] = 0.0;
       
       phi[j + i * ny][Im] = 0.0;
    }
  }
 /*
  srand(time(NULL));
  sum = 0.0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      random_num[j + i * ny] = (double) rand() / (double) RAND_MAX ;
      random_num[j + i * ny] = 2.0 * random_num[j + i * ny] - 1.0;
      random_num[j + i * ny] = random_num[j + i * ny] * noise * phi[j + i * ny][Re];
      sum += random_num[j + i * ny];
    }
  }

  mean = sum * one_by_nxny;
  printf("mean = %le\n",mean);

  // Composition profile  

  for(int i = 0; i < nx; i++){
    for (int j = 0; j < ny; j++) {
       if((double)(i - nx_half)* (double)(i - nx_half) + (double)(j - ny_half) * (double)(j - ny_half) <= R * R )
        comp[j + i * ny][Re] = 0.5 * (c_beta1 + c_beta2) + random_num[j + i * ny] - mean;
       else
          comp[j + i * ny][Re] = c_alpha;
     
        comp[j + i * ny][Im] = 0.0;
    }
  }
*/
 
   for(int i = 0; i < nx; i++){
      for(int j=0; j<ny; j++){
        if(phi[j + i * ny][Re] == 1.0)
            if( (j <= 2*nx/3 && j > (2*nx/3 - 100)) ||
                ( j> nx/3 && j < (nx/3 + 100)) )
                comp[j + i * ny][Re] = c_beta1;
            else
                comp[j + i * ny][Re] = c_beta2;
        else
           comp[j + i * ny][Re] = c_alpha;

           comp[j + i * ny][Im] = 0.0;
      }
  }

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
       dfdc[j + i * ny][Re] = comp[j + i * ny][Re];
       dfdc[j + i * ny][Im] = comp[j + i * ny][Im];
       dfdphi[j + i * ny][Re] = phi[j + i * ny][Re];
       dfdphi[j + i * ny][Im] = phi[j + i * ny][Im];
    }
  }

  sprintf(fn, "profile.in");  
  if (!(fp = fopen (fn, "w"))) {
	printf ("File:%s could not be opened \n", fn);
	exit (1);
  }

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      fprintf(fp,"%d\t%d\t%le\t%le\n",i, j, comp[j + i * ny][Re], phi[j + i * ny][Re]);
     }fprintf(fp,"\n");
  }
  fclose(fp);

}

void Read_Restart()
{
  FILE *fpread;
  char fr[100];

  sprintf (fr,"conf.%07d", initcount);
  fpread = fopen (fr, "r");
  if(fread (&comp[0], sizeof(double), 2 * nx * ny, fpread));
  if(fread (&phi[0], sizeof(double), 2 * nx * ny, fpread));
  fclose (fpread);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
       dfdc[j + i * ny][Re] = comp[j + i * ny][Re];
       dfdc[j + i * ny][Im] = comp[j + i * ny][Im];
       dfdphi[j + i * ny][Re] = phi[j + i * ny][Re];
       dfdphi[j + i * ny][Im] = phi[j + i * ny][Im];
    }
  }

}
