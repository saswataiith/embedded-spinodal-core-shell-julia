#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>

int main (void){
  int nx, ny, nz, dev, interface;

  FILE *fpread, *fpq;
  char fr[100],fw[100];
  char prefix[20];
  double total;
  char new[100];
  double *comp, *phi;

  int flag;

  fpread = fopen("InputParams", "r");
  fscanf(fpread,"%s%d", fr, &nx);
  fscanf(fpread,"%s%d", fr, &ny);
  fscanf(fpread,"%s%d", fr, &nz);
  fclose(fpread);

  int nx_half = nx/2;
  int ny_half = ny/2;
  int nz_half = nz/2;

  comp = (double *)malloc (nx * ny * nz * sizeof (double));
  phi = (double *)malloc (nx * ny * nz * sizeof (double));
  int i, j, initial, final, steps;

  printf("Initial count\n");
  scanf("%d", &initial);
  printf("Final count\n");
  scanf("%d", &final);
  printf("Steps\n");
  scanf("%d", &steps);

  printf("Enter 0:XY, 1:ZY, 2:ZX plane?\n");
  scanf("%d", &flag);

  if (flag == 0){
   for (int count = initial; count <=final; count += steps){

     sprintf (fr, "conf.%09d", count);
     printf ("%s\n", fr);
     fpread = fopen (fr, "r");
     fread (&comp[0], sizeof (double), nx * ny * nz, fpread);
     fread (&phi[0], sizeof (double), nx * ny * nz, fpread);
     fclose (fpread);
  
     sprintf (fw, "XY_view.%09d", count);
     fpq = fopen(fw,"w");
     int k = nz_half;
     for(int i = 0; i < nx; i++){
      for(int j = 0; j < ny; j++){
       fprintf(fpq, "%d\t%d\t%le\t%le\n", i, j, comp[k + nz * (j + i * ny)], phi[k + nz * (j + i * ny)]);  
      } fprintf(fpq,"\n");
     }
     fclose(fpq);
   }
  }

  if (flag == 1){
   for (int count = initial; count <=final; count += steps){

     sprintf (fr, "conf.%09d", count);
     printf ("%s\n", fr);
     fpread = fopen (fr, "r");
     fread (&comp[0], sizeof (double), nx * ny * nz, fpread);
     fread (&phi[0], sizeof (double), nx * ny * nz, fpread);
     fclose (fpread);
  
     sprintf (fw, "ZY_view.%09d", count);
     fpq = fopen(fw,"w");
     int i = nx_half;
     for(int k = 0; k < nz; k++){
      for(int j = 0; j < ny; j++){
       fprintf(fpq, "%d\t%d\t%le\t%le\n", j, k, comp[k + nz * (j + i * ny)], phi[k + nz * (j + i * ny)]);  
      } fprintf(fpq,"\n");
     }
     fclose(fpq);
   }
  }

  if (flag == 2){
   for (int count = initial; count <=final; count += steps){

     sprintf (fr, "conf.%09d", count);
     printf ("%s\n", fr);
     fpread = fopen (fr, "r");
     fread (&comp[0], sizeof (double), nx * ny * nz, fpread);
     fread (&phi[0], sizeof (double), nx * ny * nz, fpread);
     fclose (fpread);
  
     sprintf (fw, "ZX_view.%09d", count);
     fpq = fopen(fw,"w");
     int j = ny_half;
     for(int i = 0; i < nx; i++){
      for(int k = 0; k < nz; k++){
       fprintf(fpq, "%d\t%d\t%le\t%le\n", i, k, comp[k + nz * (j + i * ny)], phi[k + nz * (j + i * ny)]);  
      } fprintf(fpq,"\n");
     }
     fclose(fpq);
   }
  }
  
  free (comp);
  free (phi);
  return(0);
}

