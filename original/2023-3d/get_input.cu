#include"binary.h"

void Get_Input_Parameters (char *fnin, char *fnout)
{
  FILE *fpin, *fpcout;
  char param[100], fn[100];

  if (!(fpcout = fopen (fnout, "w"))) {
	printf ("File:%s could not be opened \n", fnout);
	exit (1);
  }
  fprintf (fpcout, "The name of this file is : %s \n", fnout);
  fprintf (fpcout, "Input is from            : %s \n", fnin);

  if (!(fpin = fopen (fnin, "r"))) {
	printf ("File: %s could not be opened \n", fnin);
	exit (1);
  }

  if(fscanf (fpin, "%s%d", param, &nx));
  if(fscanf (fpin, "%s%d", param, &ny));
  if(fscanf (fpin, "%s%d", param, &nz));
  if(fscanf (fpin, "%s%lf", param, &dx));
  if(fscanf (fpin, "%s%lf", param, &dy));
  if(fscanf (fpin, "%s%lf", param, &dz));
  if(fscanf (fpin, "%s%le",param, &dt));
  if(fscanf (fpin, "%s%d", param, &total_steps));
  if(fscanf (fpin, "%s%d", param, &print_steps));
  if(fscanf (fpin, "%s%d", param, &total_steps_1));
  if(fscanf (fpin, "%s%d", param, &print_steps_1));
  if(fscanf (fpin, "%s%d", param, &total_steps_2));
  if(fscanf (fpin, "%s%d", param, &print_steps_2));
  if(fscanf (fpin, "%s%lf", param, &R));
  if(fscanf (fpin, "%s%lf%s%lf%s%lf%s%lf", param, &A, param, &B, param, &chi, param, &P));
  if(fscanf (fpin, "%s%lf%s%lf%s%lf", param, &c_alpha, param, &c_beta1, param, &c_beta2));
  if(fscanf (fpin, "%s%lf", param, &A0));
  if(fscanf (fpin, "%s%lf", param, &c_0));
  if(fscanf (fpin, "%s%lf", param, &kappa_c));
  if(fscanf (fpin, "%s%lf", param, &kappa_phi));
  if(fscanf (fpin, "%s%lf", param, &mobility));
  if(fscanf (fpin, "%s%lf", param, &relax_coeff));
  if(fscanf (fpin, "%s%lf", param, &noise));
  if(fscanf (fpin, "%s%ld", param, &SEED));
  if(fscanf (fpin, "%s%lf", param, &Q));
  if(fscanf (fpin, "%s%le", param, &phicutoff));
  if(fscanf (fpin, "%s%d", param, &initflag));
  if(fscanf (fpin, "%s%d", param, &initcount));
  if(fscanf (fpin, "%s%le", param, &Tolerance));
  if(fscanf (fpin, "%s%d", param, &flag));
  if(fscanf (fpin, "%s%d%s%d%s%d", param, &threads_x, param, &threads_y, param, &threads_z));
  
  printf("Numerical Parameters\n");
  printf("nx = %d, ny = %d, nz = %d, dx = %0.2lf, dy = %0.2lf, dz = %0.2lf, dt = %le\n", 
		  nx, ny, nz, dx, dy, dz, dt);
  printf("Total number of simulation steps = %d\n", total_steps);
  printf("Intermediate simulation steps = %d\n", total_steps_1);
  printf("Radius = %lf\n", R);
  printf("A0 = %lf\n", A0);
  printf("c0 = %lf\n", c_0);
  
  printf("Bulk free energy coefficients A = %0.1lf, B = %0.1lf, chi = %0.2lf, P = %0.1lf\n", A, B, chi, P);
  printf("c_alpha = %0.1lf, c_beta1 = %0.1lf, c_beta2 = %0.1lf\n", c_alpha, c_beta1, c_beta2);
  printf("Kappa_c = %0.1lf, kappa_phi = %0.1lf\n", kappa_c, kappa_phi);
  printf("Relax_coeff = %le\n", relax_coeff);
  printf("mobility = %lf\n", mobility);
  printf("Noise_level = %lf, SEED = %ld\n", noise, SEED);
  printf("Tolerance = %le\n", Tolerance);
  printf("phi cut-off used = %le\n", phicutoff);
  printf("Factor for making semi-implicit = %lf\n", Q);
  printf("File read from steps = %d\n", initcount);
  printf("FFTW_FLAG = %d\n", initflag);

  if (flag == 1) printf("Homogeneous Particle\n");
  if (flag == 2) printf("Artificial Janus Particle\n");

  if(initflag == 0) {
    printf("Configuration initialized by me\n");
    initcount = 0;
  } else {
    sprintf(fn,"conf.%09d", initcount);     
    printf("Configuration read from file %s\n",fn);
  }
  //	printf("FFTW_FLAG = %d\n",fftw_flag);
  fclose (fpin);

  //fprintf (fpcout, "dev %d\n", dev);
  fprintf (fpcout, "nx %d\n", nx);
  fprintf (fpcout, "ny %d\n", ny);
  fprintf (fpcout, "nz %d\n", nz);
  fprintf (fpcout, "dx %lf\n",dx);
  fprintf (fpcout, "dy %lf\n",dy);
  fprintf (fpcout, "dz %lf\n",dz);
  fprintf (fpcout, "dt %lf\n",dt);
  fprintf (fpcout, "num_steps %d\n", total_steps);
  fprintf (fpcout, "A   %3.2f\tB  %3.2f\t  chi  %3.2f\t   P  %3.2f\n", A, B, chi, P);
  fprintf(fpcout,  "A0 = %lf\n", A0);
  fprintf (fpcout, "c_beta1 %lf\t  c_beta2  %lf\t  " 
		  "c_alpha  %lf   c_0  %lf\n",c_beta1, c_beta2, c_alpha, c_0);
  fprintf (fpcout, "kappa_c %lf\n",kappa_c);
  fprintf (fpcout, "kappa_phi %lf\n",kappa_phi);
  fprintf (fpcout, "Relax_coeff %le\n", relax_coeff);
  fprintf (fpcout, "mobility %le\n", mobility);
  fprintf (fpcout, "initcount %d\n", initcount);
  fprintf (fpcout, "initflag %d\n", initflag);
  fprintf (fpcout, "Noise_level %lf  SEED %ld\n", noise, SEED);
  fprintf (fpcout, "Semi-implicit facotr %lf\n", Q);
  fprintf (fpcout, "Tolerance %le\n", Tolerance);
  fprintf (fpcout, "Noise %le\n", noise);
  fprintf (fpcout, "phi cut-off used %le\n", phicutoff);
  if (flag == 1) fprintf(fpcout, "Homogeneous Particle, flag:%d\n", flag);
  if (flag == 2) fprintf(fpcout, "Artificial Janus Particle, flag:%d\n",flag);

  fclose(fpcout);
}
