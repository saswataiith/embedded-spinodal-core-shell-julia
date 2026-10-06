#include"binary.h"
void Output_Conf (long int steps)
{
 FILE *fpt;
  printf("%ld\n", steps); 
 char fn[100];

 sprintf (fn, "conf.%07ld", steps);

 fpt = fopen (fn, "w");
 fwrite (&dfdc[0][0], sizeof(double), 2 * nx * ny, fpt);
 fwrite (&dfdphi[0][0], sizeof(double), 2 * nx * ny , fpt);
 fclose (fpt);

/*
 sprintf (fn, "prof_gp.%09ld", steps);
 fpt = fopen (fn, "w");

 for (int i = 0; i < nx; i++) {
   for (int j = 0; j < ny; j++) {
      fprintf(fpt,"%d\t%d\t%le\t%le\n", i, j, dfdc[j + i * ny][Re], dfdphi[j + i * ny][Re]);
 }
 fprintf(fpt,"\n");
 }fclose(fpt);
*/ 

}
