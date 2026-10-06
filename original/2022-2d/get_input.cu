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

  if(fscanf (fpin, "%s%d", param, &dev));
  if(fscanf (fpin, "%s%d", param, &nx));
  if(fscanf (fpin, "%s%d", param, &ny));
  if(fscanf (fpin, "%s%lf", param, &dx));
  if(fscanf (fpin, "%s%lf", param, &dy));
  if(fscanf (fpin, "%s%le",param, &dt));
  if(fscanf (fpin, "%s%d", param, &total_steps));
  if(fscanf (fpin, "%s%d", param, &print_steps));
  if(fscanf (fpin, "%s%d", param, &total_steps_1));
  if(fscanf (fpin, "%s%d", param, &print_steps_1));
  if(fscanf (fpin, "%s%d", param, &total_steps_2));
  if(fscanf (fpin, "%s%d", param, &print_steps_2));
  if(fscanf (fpin, "%s%lf%s%lf%s%lf", param, &R, param, &Ro, param, &pre_factor));
  if(fscanf (fpin, "%s%lf%s%lf", param, &ch_ht, param, &ch_st));
  if(fscanf (fpin, "%s%lf%s%lf", param, &major, param, &minor));
  if(fscanf (fpin, "%s%lf%s%lf%s%lf%s%lf", param, &A, param, &B, param, &chi, param, &P));
  if(fscanf (fpin, "%s%lf%s%lf%s%lf", param, &c_alpha, param, &c_beta1, param, &c_beta2));
  if(fscanf (fpin, "%s%lf", param, &c_0));
  if(fscanf (fpin, "%s%lf", param, &kappa_c));
  if(fscanf (fpin, "%s%lf", param, &kappa_phi));
  if(fscanf (fpin, "%s%lf%s%lf", param, &mobility, param, &mob_tol));
  if(fscanf (fpin, "%s%lf", param, &relax_coeff));
  if(fscanf (fpin, "%s%lf", param, &noise));
  if(fscanf (fpin, "%s%d", param, &noise_steps));
  if(fscanf (fpin, "%s%d", param, &NatEveryStp));
  if(fscanf (fpin, "%s%llu", param, &SEED));
  if(fscanf (fpin, "%s%lf", param, &Q));
  if(fscanf (fpin, "%s%d", param, &flag));
  if(fscanf (fpin, "%s%le", param, &phicutoff));
  if(fscanf (fpin, "%s%d", param, &initflag));
  if(fscanf (fpin, "%s%d", param, &initcount));
  if(fscanf (fpin, "%s%le", param, &Tolerance));
  if(fscanf (fpin, "%s%d%s%d", param, &threads_x, param, &threads_y));
  /*if(fscanf (fpin, "%s%d", param, &fftw_flag));
  */
  printf("Numerical Parameters\n");
  printf("GPU used = %d\n", dev);
  printf("nx = %d, ny = %d, dx = %0.2lf, dy = %0.2lf, dt = %le\n", 
		  nx, ny, dx, dy, dt);
  printf("Total number of simulation steps = %d\n", total_steps);
  printf("Intermediate simulation steps = %d\n", total_steps_1);
  printf("Radius = %lf, phi width prefactor = %lf\n", R, pre_factor);
  if(flag == 0){ printf ("Circle Radius: %0.3lf\n", R);}
  if(flag == 1){ printf ("Channel: length: %d, width: %0.3lf\n", nx, ch_ht);}
  if(flag == 2){ printf ("Ellipse: major axis %0.3lf, minor axis %0.3lf\n", major, minor);}
  if(flag == 3){ printf ("Circle Inner Radius: %0.3lf, Outer Radius = %0.3lf\n", R, Ro);}
  if(flag == 4){ printf ("Circle inside Ellipse, Cr Radius: %0.3lf, major axis %0.3lf, minor axis %0.3f, Eq radius%0.3lf\n", R, major, minor, Ro);}
  if(flag == 5){ printf ("Ellipse inside Cirlce, Cr Radius: %0.3lf, major axis %0.3lf, minor axis %0.3f, Eq radius%0.3lf\n", Ro, major, minor, R);}
  printf("Bulk free energy coefficients A = %0.1lf, B = %0.1lf, chi = %lf, P = %0.1lf\n", 
		  A,B, chi, P);
  printf("c_alpha = %0.1lf, c_beta1 = %0.1lf, c_beta2 = %0.1lf\n", 
		  c_alpha, c_beta1, c_beta2);
  printf("Kappa_c = %0.1lf, kappa_phi = %0.1lf\n", kappa_c, kappa_phi);
  printf("Mobility = %0.1lf, mob_tol = %lf, relax_coeff = %le\n", mobility, mob_tol, relax_coeff);
  printf("Noise_level = %lf, SEED = %llu\n", noise, SEED);
  printf("Tolerance = %le\n", Tolerance);
  printf("Noise added for %d steps and later at every %d step\n", noise_steps, NatEveryStp);
  printf("phi cut-off used = %le\n", phicutoff);
  printf("Factor for making semi-implicit = %lf\n", Q);
  printf("File read from steps = %d\n", initcount);
  printf("FFTW_FLAG = %d\n", initflag);
  if(initflag == 0) {
    printf("Configuration initialized by me\n");
    initcount = 0;
  } else {
    sprintf(fn,"conf.%11d", initcount);     
    printf("Configuration read from file %s\n",fn);
  }
  //	printf("FFTW_FLAG = %d\n",fftw_flag);
  fclose (fpin);

  //fprintf (fpcout, "dev %d\n", dev);
  fprintf (fpcout, "nx %d\n", nx);
  fprintf (fpcout, "ny %d\n", ny);
  fprintf (fpcout, "dx %lf\n",dx);
  fprintf (fpcout, "dy %lf\n",dy);
  fprintf (fpcout, "dt %lf\n",dt);
  fprintf (fpcout, "num_steps %d\n", total_steps);
  fprintf (fpcout, "A   %3.2f\tB  %3.2f\t  chi  %3.2f\t   P  %3.2f\n", A, B, chi, P);
  fprintf (fpcout, "c_beta1 %lf\t  c_beta2  %lf\t  " 
		 "c_alpha  %lf   c_0  %lf\n",c_beta1, c_beta2, c_alpha, c_0);
  fprintf (fpcout, "kappa_c %lf\n",kappa_c);
  fprintf (fpcout, "kappa_phi %lf\n",kappa_phi);
  fprintf (fpcout, "mobility %lf\n", mobility);
  fprintf (fpcout, "relax_coeff %lf\n", relax_coeff);
  fprintf (fpcout, "initcount %d\n", initcount);
  fprintf (fpcout, "initflag %d\n", initflag);
  fprintf (fpcout, "Noise_level %lf  SEED %llu\n", noise, SEED);
  fprintf (fpcout, "Semi-implicit facotr %lf\n", Q);
  fprintf (fpcout, "Tolerance %le\n", Tolerance);
  fprintf (fpcout, "Shape flag %d\n", flag);
  fprintf (fpcout, "Noise %le, noise_steps %d\n", noise, noise_steps);
  fprintf (fpcout, "phi cut-off used %le\n", phicutoff);
  if(flag == 0){ fprintf (fpcout, "Circle Radius: %0.3lf\n", R);}
  if(flag == 1){ fprintf (fpcout, "Channel: length: %d, width: %0.3lf\n", nx, ch_ht);}
  if(flag == 2){ fprintf (fpcout, "Ellipse: major axis %0.3lf, minor axis %0.3lf\n", major, minor);}
  if(flag == 3){ fprintf (fpcout, "Circle Inner Radius: %0.3lf, Outer Radius: %0.3lf\n", R, Ro);}
  if(flag == 4){ fprintf (fpcout, "Circle inside Ellipse, Cr Radius: %0.3lf, major axis %0.3lf, minor axis %0.3f, Eq radius%0.3lf\n", R, major, minor, Ro);}
  if(flag == 5){ fprintf (fpcout, "Ellipse inside Cirlce, Cr Radius: %0.3lf, major axis %0.3lf, minor axis %0.3f, Eq radius%0.3lf\n", Ro, major, minor, R);}
  fclose(fpcout);
}
