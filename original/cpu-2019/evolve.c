#include"binary.h"
void Evolve()//, fftw_complex *dfdc, fftw_complex *dfdphi)
{
   void Output_Conf(long int steps);
   int loop_condition, count;
   
   double dkx, kx, dky, ky, kpow2, kpow4;

   double rc, fp, rc_new, rphi_new;
   double rphi, fpphi;

   double total, err, f, sum, mean;
   double *temprealA, *temprealB;
   double lhs, rhs, lhse, rhse;

   double errA, maxerrorA;
   double errB, maxerrorB;

   temprealA = (double *) malloc(sizeof(double) * nx * ny);
   temprealB = (double *) malloc(sizeof(double) * nx * ny);

   dkx = 2.0 * PI / ((double) nx * dx);
   dky = 2.0 * PI / ((double) ny * dy);
  // for (int i = 0; i < nx; i++) 
  //	tempreal[i] = comp[i][Re];

   loop_condition = 1;

   fftw_execute_dft(p_up, comp, comp);
   fftw_execute_dft(p_up, phi, phi);

   alloycomp = comp[0][Re] * one_by_nxny;

   clock_t begin = clock();

   //Start the evolution

   for (count = 0; count <= num_steps; count++) {
   if (((count % print_steps) == 0) || (count == num_steps)
		|| (loop_condition == 0)) {
		printf("time steps %d\n", count);
	printf("total_time=%lf\n", sim_time);
	printf("writing configuration to file!\n");

	Output_Conf(count);
     }
     if (count > num_steps || loop_condition == 0)
	 break;

     // Evaluate dfdc in real space

     double ctemp;
     double ptemp;

     double gphi, hphi, hprime, gprime;
     for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
	     ctemp = dfdc[j + i * ny][Re];	//this dfdc is composition
	     ptemp = dfdphi[j + i * ny][Re];	//this dfdphi is phi
	     hphi = ptemp * ptemp * ptemp * (10.0 - 15.0 * ptemp + 6.0 * ptemp * ptemp);
	     hprime = 30.0 * (ptemp * ptemp - 2.0 * ptemp * ptemp * ptemp + ptemp * ptemp * ptemp * ptemp);
	     gphi = (ptemp * ptemp) * (1.0 - ptemp) * (1.0 - ptemp);
	     gprime = 2.0 * ptemp - 6.0 * ptemp * ptemp + 4.0 * ptemp * ptemp * ptemp;

	     dfdc[j + i * ny][Re] = 2.0 * A * (1.0 - hphi) * (ctemp - c_alpha) 
                          + 2.0 * B * hphi * (ctemp - c_beta1) * (ctemp - c_beta2) 
			  * (2.0 *ctemp - c_beta1 - c_beta2) - chi * P * gphi ;
	     dfdc[j + i * ny][Im] = 0.0;	//this dfdc is derivative of f wrt c

	     dfdphi[j + i * ny][Re] = -1.0 * hprime * A * (ctemp - c_alpha) * (ctemp - c_alpha) 
		             + hprime * B * (ctemp - c_beta1) * (ctemp - c_beta1) 
			      * (ctemp - c_beta2) * (ctemp - c_beta2) + (1.0 - chi * ctemp ) * P * gprime ;
	     dfdphi[j + i * ny][Im] = 0.0;	// this dfdphi is derivative of f wrt phi
        }
     }
FILE *po;
po = fopen("check","w");
     for (int i = 0; i < nx; i++) {
       for (int j = 0; j < ny; j++) {
            fprintf(po,"%d\t%d\t%le\t%le\n",i, j, dfdc[j + i * ny][Re], dfdphi[j + i * ny][Re]);
       } fprintf(po,"\n");
     }fclose(po);

     fftw_execute_dft(p_up, dfdc, dfdc);
     fftw_execute_dft(p_up, dfdphi, dfdphi);

     for (int i = 0; i < nx; i++) {
       for (int j = 0; j < ny; j++) {
  	  if (i <= nx_half)
	     kx = (double) i *dkx;
	  else
	     kx = (double) (i - nx) * dkx;

  	  if (j <= ny_half)
	     ky = (double) j *dky;
	  else
	     ky = (double) (j - ny) * dky;
	kpow2 = kx * kx + ky * ky;
	kpow4 = kpow2 * kpow2;

	lhs = 1.0 + 2.0 * mobility * kappa_c * kpow4 * dt;

	rc = comp[j + i * ny][Re];
	fp = dfdc[j + i * ny][Re];
   	muA[j + i * ny][Re] = - c_beta2 * (fp + 2.0 * kappa_c * kpow2 * rc);
	muB[j + i * ny][Re] = (fp + 2.0 * kappa_c * kpow2 * rc) + muA[j + i * ny][Re];
	rhs = rc - mobility * kpow2 * dt * fp;
	rc_new = rhs / lhs;
	comp[j + i * ny][Re] = rc_new;
	dfdc[j + i * ny][Re] = comp[j + i * ny][Re];

	rc = comp[j + i * ny][Im];
	fp = dfdc[j + i * ny][Im];
   	muA[j + i * ny][Im] = - c_beta2 * (fp + 2.0 * kappa_c * kpow2 * rc);
	muB[j + i * ny][Im] = (fp + 2.0 * kappa_c * kpow2 * rc) + muA[j + i * ny][Im];
	rhs = rc - mobility * kpow2 * dt * fp;
	rc_new = rhs / lhs;
	comp[j + i * ny][Im] = rc_new;
	dfdc[j + i * ny][Im] = comp[j + i * ny][Im];

	lhse = 1.0 + 2.0 * relax_coeff * kappa_phi * kpow2 * dt;

	rphi = phi[j + i * ny][Re];
	fpphi = dfdphi[j + i * ny][Re];
	rhse = rphi - relax_coeff * dt * fpphi;
	rphi_new = rhse / lhse;
	phi[j + i * ny][Re] = rphi_new;
	dfdphi[j + i * ny][Re] = phi[j + i * ny][Re];

	//printf("%lf\n", dfdphi[i][Re]);
	rphi = phi[j + i * ny][Im];
	fpphi = dfdphi[j + i * ny][Im];
	rhse = rphi - relax_coeff * dt * fpphi;
	rphi_new = rhse / lhse;
	phi[j + i * ny][Im] = rphi_new;
	dfdphi[j + i * ny][Im] = phi[j + i * ny][Im];
       }
     }

     //Check for conservation of mass
     total = dfdc[0][Re] * one_by_nxny;
     err = fabs(total - alloycomp);
     if (err > COMPERR) {
	printf("ELEMENTS ARE NOT CONSERVED,SORRY!!!!\n");
	printf("error=%lf\n", err);
	exit(0);
     }

     fftw_execute_dft(p_dn, dfdc, dfdc);
     fftw_execute_dft(p_dn, dfdphi, dfdphi);
     fftw_execute_dft(p_dn, muA, muA);
     fftw_execute_dft(p_dn, muB, muB);
     
     for (int i = 0; i < nx; i++) {
       for (int j = 0; j < ny; j++) {
	 dfdc[j + i * ny][Re] *= one_by_nxny;
	 dfdc[j + i * ny][Im] *= one_by_nxny;
	 dfdphi[j + i * ny][Re] *= one_by_nxny;
	 dfdphi[j + i * ny][Im] *= one_by_nxny;
	 muA[j + i * ny][Re] *= one_by_nxny;
	 muA[j + i * ny][Im] *= one_by_nxny;
	 muB[j + i * ny][Re] *= one_by_nxny;
	 muB[j + i * ny][Im] *= one_by_nxny;
       }
     }

     //Check for bounds
     for (int i = 0; i < nx; i++) {
       for (int j = 0; j < ny; j++) {
	 if (dfdc[j + i * ny][Re] < -0.2 || dfdc[j + i * ny][Re] > 1.2) {
	    printf("Compositions out of bounds. Exiting\n");
	    exit(0);
	 }
       }
     }
     for (int i = 0; i < nx; i++) {
       for (int j = 0; j < ny; j++) {
         f = A * (dfdc[j + i * ny][Re] - c_alpha) * (dfdc[j + i * ny][Re] - c_alpha) * (1.0 - hphi) 
                       + (B * (dfdc[j + i * ny][Re] - c_beta1) * (dfdc[j + i * ny][Re] - c_beta1) 
		       * (dfdc[j + i * ny][Re] - c_beta2) * (dfdc[j + i * ny][Re] - c_beta2 )) * hphi
	 	       + (1.0 - chi * dfdc[j + i * ny][Re]) * P * dfdphi[j + i * ny][Re] * dfdphi[j + i * ny][Re] 
		       * (1.0 - dfdphi[j + i * ny][Re]) * (1.0 - dfdphi[j + i * ny][Re]);
         muA[j + i * ny][Re] += f;
         muB[j + i * ny][Re] += f;
       }
     }

     if (count > 0 ){   
       //Check for convergence
       maxerrorA = 0.0;
       maxerrorB = 0.0;
       for (int i = 0; i < nx; i++) {
         for (int j = 0; j < ny; j++) {
 	    errA = fabs(temprealA[j + i * ny] - muA[j + i * ny][Re]);
	    errB = fabs(temprealB[j + i * ny] - muB[j + i * ny][Re]);
	    if (errA > maxerrorA)
	       maxerrorA = errA;
	    if (errB > maxerrorB)
	       maxerrorB = errB;
         }
       }

       if ( (maxerrorA <= Tolerance) || (maxerrorB <= Tolerance) ) {
  	  printf("maxerrorA=%lf\tmaxerrorB=%lf\tnumbersteps=%d\n", maxerrorA, maxerrorA, count);
          loop_condition = 0;
       }
     }
     sim_time = sim_time + dt;
     for (int i = 0; i < nx; i++) {
         for (int j = 0; j < ny; j++) {
	   temprealA[j + i * ny] = muA[j + i * ny][Re]; 
	   temprealB[j + i * ny] = muB[j + i * ny][Re]; 
         }
     }
     /*
     sum = 0.0; 
     if (count % print_steps == 0){
        for (int i = 0; i < nx; i++)
            sum += dfdphi[i][Re]; 
        
        mean = sum * one_by_nxny;
        printf("meanPhi = %le\n", mean);
     }
     */
   }
   clock_t end = clock();
   double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
   printf("Time consumed in 'count' loop = %lf\n", time_spent);
}
