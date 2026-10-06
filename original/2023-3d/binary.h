
#ifndef BINARY_H
#define BINARY_H

#include<stdio.h>
#include<stdlib.h>
#include<cuda.h>
#include<curand.h>
#include"check.h"
#include<sys/time.h>
#include</usr/local/cuda-11.2/targets/x86_64-linux/include/cufft.h>
#include"cub/cub.cuh"

#define COMPERR 1.0e-08
double Tolerance;

double *d_gradmu;
int flag, horizontal, vertical, angle;
int dev;
long SEED;
int counter;
/*
 * Free energy(thermodynamics) parameters
 */
   cufftDoubleComplex *phi, *comp;

   cufftDoubleComplex *d_phi, *d_comp;
   cufftDoubleComplex *d_dfdc, *d_dfdphi;

/*
 * Memory Parameters
 */
   size_t nBytes;
   size_t nBytesComplex;

/*
 * Configuration to be initialized or to be read
 */
   int initcount, initflag;

  double mobility;
/*
 * Geometeric Parameters
 */
   int nx, ny, nz;
   int nx_half, ny_half, nz_half;
   double one_by_nxnynz;
   double one_by_counter;
   double dx, dy, dz;
   double R, Ro, Ri;
   double noise, SusNoise, mean_noise, sd;
   double Q; //Factor to make discretization semi-implicit
/* 
 * Random Number Generator Parameters
 */
   double *d_randm, *randm;
   curandGenerator_t gen;
/*
 * kx and ky
 */
   double *kx1, *kx2, *ky1, *ky2, *kz1, *kz2;
   double *d_kx1, *d_kx2, *d_ky1, *d_ky2, *d_kz1, *d_kz2;

/*
 * Free energy Parameters
*/
   double A, B, chi, P, c_0, A0;
   double c_alpha, c_beta1, c_beta2;

/*
 * Kinetic Parameters
 */ 
   double kappa_c, kappa_phi;
   double mp, mm, ms, relax_coeff;

/*
 * Time Paramters
 */ double dt;
    double sim_time, total_time;

/*
 * Printing parameters
 */ //Total Number of simulation steps
    int total_steps, total_steps_1, total_steps_2;
    //step interval for printing output
    int print_steps, print_steps_1, print_steps_2;
    //noise steps
    int noise_steps, noise_Atevery_steps;

    double pre_factor, sum, mean;

    double phicutoff;

    int Nofpart;

    FILE *fpout;

int threads_x, threads_y, threads_z;
int blocks_x, blocks_y, blocks_z;

#endif
