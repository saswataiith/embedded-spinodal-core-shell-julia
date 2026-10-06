__global__ void kernel_cphi_der(cufftDoubleComplex *d_dfdc, 
		cufftDoubleComplex *d_dfdphi, double mobility,
		double A, double A0, double B, double chi, double P, 
		double c_alpha, double c_beta1, double c_beta2, 
		int nx, int ny, int nz){

  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
  unsigned int tid = iz + nz * (iy + ix * ny);
  
  double  ctemp, ptemp, hphi, gphi, hprime, gprime;
  
  if(tid < nx * ny * nz){

     ctemp = d_dfdc[tid].x;
     ptemp = d_dfdphi[tid].x;
     hphi = ptemp * ptemp * ptemp * (10.0 - 15.0 * ptemp + 6.0 * ptemp * ptemp);
     gphi = (ptemp * ptemp ) * (1.0 - ptemp) * (1.0 - ptemp);
     hprime = 30.0 * (ptemp * ptemp - 2.0 * ptemp * ptemp * ptemp +
		  ptemp * ptemp * ptemp * ptemp); 
     gprime = 2.0 * ptemp - 6.0 * ptemp *ptemp + 4.0 * ptemp * ptemp * ptemp;

     d_dfdc[tid].x = 2.0 * A * (1.0 - hphi) * (ctemp - c_alpha) + 
	             2.0 * B * hphi * (ctemp - c_beta1) * (ctemp - c_beta2) *
		     (2.0 * ctemp - c_beta1 - c_beta2) - chi * P * gphi;
     d_dfdc[tid].y = 0.0;
     
     d_dfdphi[tid].x = -1.0 * hprime * (A * (ctemp - c_alpha) * (ctemp - c_alpha) + A0)
                       + hprime * B * (ctemp - c_beta1) * (ctemp - c_beta1) * 
                	 (ctemp - c_beta2) * (ctemp - c_beta2) + (1.0 - chi * ctemp) * P * gprime;

     d_dfdphi[tid].y = hphi * mobility;

  }
}
    
__global__ void Compute_mu(cufftDoubleComplex *d_muz, cufftDoubleComplex *d_comp, 
		           cufftDoubleComplex *d_dfdc,
			 double *d_kx2, double *d_ky2, double *d_kz2,
				double kappa_c, int nx, int ny, int nz){
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
  unsigned int tid = iz + nz * (iy + ix * ny);
 
  if(tid < nx * ny * nz){
    d_muz[tid].x = d_dfdc[tid].x + 
	    	     2.0 * kappa_c * (d_kx2[ix] * d_kx2[ix] + d_ky2[iy] * d_ky2[iy] + d_kz2[iz] * d_kz2[iz]) * d_comp[tid].x;

    d_muz[tid].y = d_dfdc[tid].y + 
	             2.0 * kappa_c * (d_kx2[ix] * d_kx2[ix] + d_ky2[iy] * d_ky2[iy] + d_kz2[iz] * d_kz2[iz]) * d_comp[tid].y;

    d_dfdc[tid].x = 0.0;
    d_dfdc[tid].y = 0.0;
  }

}

__global__ void Muscaling(cufftDoubleComplex *d_mu, double one_by_nxnynz, 
				int nx, int ny, int nz){

  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
  unsigned int tid = iz + nz * (iy + ix * ny);
     
  if(tid < nx * ny * nz){
     d_mu[tid].x *= one_by_nxnynz;
     d_mu[tid].y = 0.0;
  }
}

__global__ void ComputeGradmu_x(double *d_mux,
				cufftDoubleComplex *d_mu,
				double one_by_nxnynz,
				int nx, int ny, int nz, double dx)
{
    unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
    unsigned int k = threadIdx.z + blockIdx.z * blockDim.z;

    unsigned int idx = k + nz * (j + i * ny);

    int xp[2], xm[2];

    if (i < nx && j < ny && k < nz)
    {
        if (i == 0)
        {
            xp[0] = (j + (i + 1) * ny) * nz + k;
            xp[1] = (j + (i + 2) * ny) * nz + k;
            xm[0] = (j + (nx - 1) * ny) * nz + k;
            xm[1] = (j + (nx - 2) * ny) * nz + k;
        }
        else if (i == 1)
        {
            xp[0] = (j + (i + 1) * ny) * nz + k;
            xp[1] = (j + (i + 2) * ny) * nz + k;
            xm[0] = j * nz + k;
            xm[1] = (j + (nx - 1) * ny) * nz + k;
        }
        else if (i == nx - 2)
        {
            xp[0] = (j + (i + 1) * ny) * nz + k;
            xp[1] = j * nz + k;
            xm[0] = (j + (i - 1) * ny) * nz + k;
            xm[1] = (j + (i - 2) * ny) * nz + k;
        }
        else if (i == nx - 1)
        {
            xp[0] = j * nz + k;
            xp[1] = (j + ny) * nz + k;
            xm[0] = (j + (i - 1) * ny) * nz + k;
            xm[1] = (j + (i - 2) * ny) * nz + k;
        }
        else
        {
            xp[0] = (j + (i + 1) * ny) * nz + k;
            xp[1] = (j + (i + 2) * ny) * nz + k;
            xm[0] = (j + (i - 1) * ny) * nz + k;
            xm[1] = (j + (i - 2) * ny) * nz + k;
        }

        d_mux[idx] = (-1 * d_mu[xp[1]].x + 8 * d_mu[xp[0]].x - 8 * d_mu[xm[0]].x + d_mu[xm[1]].x) / (12 * dx);
        //d_mux[idx] = (d_mu[xp[0]].x - d_mu[xm[0]].x) / (2.0 * dx);

    }
    __syncthreads();
}

__global__ void ComputeGradmu_y(double *d_muy,
				cufftDoubleComplex *d_mu,
				double one_by_nxnynz,
				int nx, int ny, int nz, double dy)
{
    unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
    unsigned int k = threadIdx.z + blockIdx.z * blockDim.z;

    int yp[2], ym[2];

    unsigned int idx = k + nz * (j + i * ny);

    if (i < nx && j < ny && k <  nz)
    {
        if (j == 0)
        {
            yp[0] = (j + 1 + i * ny) * nz + k;
            yp[1] = (j + 2 + i * ny) * nz + k;
            ym[0] = (ny - 1 + i * ny) * nz + k;
            ym[1] = (ny - 2 + i * ny) * nz + k;
        }
        else if (j == 1)
        {
            yp[0] = (j + 1 + i * ny) * nz + k;
            yp[1] = (j + 2 + i * ny) * nz + k;
            ym[0] = i * ny * nz + k;
            ym[1] = (ny - 1 + i * ny) * nz + k;
        }
        else if (j == ny - 2)
        {
            yp[0] = (j + 1 + i * ny) * nz + k;
            yp[1] = i * ny * nz + k;
            ym[0] = (j - 1 + i * ny) * nz + k;
            ym[1] = (j - 2 + i * ny) * nz + k;
        }
        else if (j == ny - 1)
        {
            yp[0] = i * ny * nz + k;
            yp[1] = (1 + i * ny) * nz + k;
            ym[0] = (j-1 + i * ny) * nz + k;
            ym[1] = (j-2 + i * ny) * nz + k;
        }
        else
        {
            yp[0] = (j + 1 + i * ny) * nz + k;
            yp[1] = (j + 2 + i * ny) * nz + k;
            ym[0] = (j - 1 + i * ny) * nz + k;
            ym[1] = (j - 2 + i * ny) * nz + k;
        }

        d_muy[idx] = (-1 * d_mu[yp[1]].x + 8 * d_mu[yp[0]].x - 8 * d_mu[ym[0]].x + d_mu[ym[1]].x)/(12 * dy);
        //d_muy[idx] = (d_mu[yp[0]].x - d_mu[ym[0]].x) / (2.0 * dy); 
    }
    __syncthreads();
}

__global__ void ComputeGradmu_z(double *d_muz,
				cufftDoubleComplex *d_mu,
				double one_by_nxnynz,
				int nx, int ny, int nz, double dz)
{
    unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
    unsigned int k = threadIdx.z + blockIdx.z * blockDim.z;

    int zp[2], zm[2];

    unsigned int idx = k + nz * (j + i * ny);

    if (i < nx && j < ny && k <  nz)
    {
        if (k == 0)
        {
            zp[0] = (j + i * ny ) * nz + k+1;
            zp[1] = (j + i * ny ) * nz + k+2;
            zm[0] = (j + i * ny ) * nz +  nz-1;
            zm[1] = (j + i * ny ) * nz +  nz-2;
        }
        else if (k == 1)
        {
            zp[0] = (j + i * ny ) * nz + k+1;
            zp[1] = (j + i * ny ) * nz + k+2;
            zm[0] = (j + i * ny ) * nz;
            zm[1] = (j + i * ny ) * nz +  nz-1;
        }
        else if (k == nz - 2)
        {
            zp[0] = (j + i * ny ) * nz + k+1;
            zp[1] = (j + i * ny ) * nz;
            zm[0] = (j + i * ny ) * nz + k-1;
            zm[1] = (j + i * ny ) * nz + k-2;
        }
        else if (k == nz - 1)
        {
            zp[0] = (j + i * ny ) * nz;
            zp[1] = (j + i * ny ) * nz + 1;
            zm[0] = (j + i * ny ) * nz + k-1;
            zm[1] = (j + i * ny ) * nz + k-2;
        }
        else
        {
            zp[0] = (j + i * ny ) * nz + k+1;
            zp[1] = (j + i * ny ) * nz + k+2;
            zm[0] = (j + i * ny ) * nz + k-1;
            zm[1] = (j + i * ny ) * nz + k-2;
        }

        d_muz[idx] = (-1 * d_mu[zp[1]].x + 8*d_mu[zp[0]].x - 8*d_mu[zm[0]].x + d_mu[zm[1]].x)/(12 * dz);
        //d_muz[idx] = (d_mu[zp[0]].x - d_mu[zm[0]].x) / (2.0 * dz); 
    }
    __syncthreads();
}

__global__ void ComputeGrad_x(cufftDoubleComplex *d_dfdc,
		              double *d_mux, cufftDoubleComplex *d_dfdphi, 
			   	int nx, int ny, int nz, double dx)
{
    unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
    unsigned int k = threadIdx.z + blockIdx.z * blockDim.z;

    unsigned int idx = (j + i *  ny )* nz + k;

    int xp[2], xm[2];
 
    if (i < nx && j <  ny  && k <  nz)
    {
        if (i == 0)
        {
            xp[0] = (j + (i+1)* ny )* nz + k;
            xp[1] = (j + (i+2)* ny )* nz + k;
            xm[0] = (j + (nx-1)* ny )* nz + k;
            xm[1] = (j + (nx-2)* ny )* nz + k;
        }
        else if (i == 1)
        {
            xp[0] = (j + (i+1)* ny )* nz + k;
            xp[1] = (j + (i+2)* ny )* nz + k;
            xm[0] = j* nz + k;
            xm[1] = (j + (nx-1)* ny )* nz + k;
        }
        else if (i == nx - 2)
        {
            xp[0] = (j + (i+1)* ny )* nz + k;
            xp[1] = j* nz + k;
            xm[0] = (j + (i-1)* ny )* nz + k;
            xm[1] = (j + (i-2)* ny )* nz + k;
        }
        else if (i == nx - 1)
        {
            xp[0] = j* nz + k;
            xp[1] = (j +  ny )* nz + k;
            xm[0] = (j + (i-1)* ny )* nz + k;
            xm[1] = (j + (i-2)* ny )* nz + k;
        }
        else
        {
            xp[0] = (j + (i+1)* ny )* nz + k;
            xp[1] = (j + (i+2)* ny )* nz + k;
            xm[0] = (j + (i-1)* ny )* nz + k;
            xm[1] = (j + (i-2)* ny )* nz + k;
        }

	d_dfdc[idx].x += (-1.0 * d_dfdphi[xp[1]].y * d_mux[xp[1]] + 8.0 * d_dfdphi[xp[0]].y * d_mux[xp[0]]
			  -8.0 * d_dfdphi[xm[0]].y * d_mux[xm[0]] + d_dfdphi[xm[1]].y * d_mux[xm[1]]) / (12.0 * dx);

	//d_dfdc[idx].x += (d_dfdphi[xp[0]].y * d_mux[xp[0]] - d_dfdphi[xm[0]].y * d_mux[xm[0]]) / (2.0 * dx);
    }
    __syncthreads();
}

__global__ void ComputeGrad_y(cufftDoubleComplex *d_dfdc,
				double *d_muy, cufftDoubleComplex *d_dfdphi, 
	       			int nx, int ny, int nz, double dy)
{
    unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
    unsigned int k = threadIdx.z + blockIdx.z * blockDim.z;

    unsigned int idx = (j + i * ny) * nz + k;

    int yp[2], ym[2];

    if (i < nx && j < ny && k < nz)
    {
        if (j == 0)
        {
            yp[0] = (j + 1 + i * ny ) * nz + k;
            yp[1] = (j + 2 + i * ny ) * nz + k;
            ym[0] = ( ny - 1 + i * ny ) * nz + k;
            ym[1] = ( ny - 2 + i * ny ) * nz + k;
        }
        else if (j == 1)
        {
            yp[0] = (j+1 + i* ny )* nz + k;
            yp[1] = (j+2 + i* ny )* nz + k;
            ym[0] = i* ny * nz + k;
            ym[1] = ( ny -1 + i* ny )* nz + k;
        }
        else if (j ==  ny  - 2)
        {
            yp[0] = (j+1 + i* ny )* nz + k;
            yp[1] = i* ny * nz + k;
            ym[0] = (j-1 + i* ny )* nz + k;
            ym[1] = (j-2 + i* ny )* nz + k;
        }
        else if (j ==  ny  - 1)
        {
            yp[0] = i* ny * nz + k;
            yp[1] = (1 + i* ny )* nz + k;
            ym[0] = (j-1 + i* ny )* nz + k;
            ym[1] = (j-2 + i* ny )* nz + k;
        }
        else
        {
            yp[0] = (j+1 + i* ny )* nz + k;
            yp[1] = (j+2 + i* ny )* nz + k;
            ym[0] = (j-1 + i* ny )* nz + k;
            ym[1] = (j-2 + i* ny )* nz + k;
        }

	d_dfdc[idx].x += (-1.0 * d_dfdphi[yp[1]].y * d_muy[yp[1]] + 8.0 * d_dfdphi[yp[0]].y * d_muy[yp[0]]
 			  -8.0 * d_dfdphi[ym[0]].y * d_muy[ym[0]] + d_dfdphi[ym[1]].y * d_muy[ym[1]]) / (12.0 * dy);

	//d_dfdc[idx].x += (d_dfdphi[yp[0]].y * d_muy[yp[0]] - d_dfdphi[ym[0]].y * d_muy[ym[0]]) / (2.0 * dy);

    }
    __syncthreads();
}

__global__ void ComputeGrad_z(cufftDoubleComplex *d_dfdc,
				double *d_muz, cufftDoubleComplex *d_dfdphi, 
	       			int nx, int ny, int nz, double dz)
{
    unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
    unsigned int k = threadIdx.z + blockIdx.z * blockDim.z;

    unsigned int idx = (j + i *  ny )* nz + k;

    int zp[2], zm[2];

    if (i < nx && j <  ny  && k <  nz)
    {
        if (k == 0)
        {
            zp[0] = (j + i* ny )* nz + k+1;
            zp[1] = (j + i* ny )* nz + k+2;
            zm[0] = (j + i* ny )* nz +  nz-1;
            zm[1] = (j + i* ny )* nz +  nz-2;
        }
        else if (k == 1)
        {
            zp[0] = (j + i* ny )* nz + k+1;
            zp[1] = (j + i* ny )* nz + k+2;
            zm[0] = (j + i* ny )* nz;
            zm[1] = (j + i* ny )* nz +  nz-1;
        }
        else if (k ==  nz - 2)
        {
            zp[0] = (j + i* ny )* nz + k+1;
            zp[1] = (j + i* ny )* nz;
            zm[0] = (j + i* ny )* nz + k-1;
            zm[1] = (j + i* ny )* nz + k-2;
        }
        else if (k ==  nz - 1)
        {
            zp[0] = (j + i* ny )* nz;
            zp[1] = (j + i* ny )* nz + 1;
            zm[0] = (j + i* ny )* nz + k-1;
            zm[1] = (j + i* ny )* nz + k-2;
        }
        else
        {
            zp[0] = (j + i* ny )* nz + k+1;
            zp[1] = (j + i* ny )* nz + k+2;
            zm[0] = (j + i* ny )* nz + k-1;
            zm[1] = (j + i* ny )* nz + k-2;
        }

	d_dfdc[idx].x += (-1.0 * d_dfdphi[zp[1]].y * d_muz[zp[1]] + 8.0 * d_dfdphi[zp[0]].y * d_muz[zp[0]]
        			  -8.0 * d_dfdphi[zm[0]].y * d_muz[zm[0]] + d_dfdphi[zm[1]].y * d_muz[zm[1]]) / (12.0 * dz);
       
	//d_dfdc[idx].x += (d_dfdphi[zp[0]].y * d_xyz[zp[0]] - d_dfdphi[zm[0]].y * d_xyz[zm[0]]) / (2.0 * dz);
    }
    __syncthreads();
}

__global__ void reset_dfdphi(cufftDoubleComplex *d_dfdphi, int nx, 
				int ny, int nz){
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
  unsigned int tid = iz + nz * (iy + ix * ny);

  if(tid < nx * ny * nz){
    d_dfdphi[tid].y = 0.0;
  }
}
__global__ void kernel_comp_phi(cufftDoubleComplex *d_comp,
	       	cufftDoubleComplex *d_dfdc, 
		cufftDoubleComplex *d_phi, cufftDoubleComplex *d_dfdphi,
       	        double *d_kx1, double *d_ky1, double *d_kz1, 
		double *d_kx2, double *d_ky2, double *d_kz2,
		double dt, double kappa_c, double kappa_phi, double Q, 
		double relax_coeff, int nx, int ny, int nz){

  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
  unsigned int tid = iz + nz * (iy + ix * ny);

  double lhs, lhse, kpow2;
 
  if(tid < nx * ny * nz){
    kpow2 = d_kx2[ix] * d_kx2[ix] + d_ky2[iy] * d_ky2[iy] + d_kz2[iz] * d_kz2[iz];
    lhs = 1.0 + 2.0 * Q * kappa_c * kpow2 * kpow2 * dt;

    d_comp[tid].x = (dt * d_dfdc[tid].x +
		     (1.0 + 2.0 * Q * kappa_c * kpow2 * kpow2 * dt) * 
		      d_comp[tid].x) / lhs;
    d_dfdc[tid].x = d_comp[tid].x;

    d_comp[tid].y = (dt * d_dfdc[tid].y +
		     (1.0 + 2.0 * Q * kappa_c * kpow2 * kpow2 * dt) * 
		      d_comp[tid].y) / lhs;
    d_dfdc[tid].y = d_comp[tid].y;
     
    lhse = 1.0 + 2.0 * relax_coeff * kappa_phi * kpow2 * dt;

    d_phi[tid].x = (d_phi[tid].x - relax_coeff * dt * d_dfdphi[tid].x) / lhse;
    d_dfdphi[tid].x = d_phi[tid].x;

    //Imag part
    d_phi[tid].y = (d_phi[tid].y - relax_coeff * dt * d_dfdphi[tid].y) / lhse;
    d_dfdphi[tid].y = d_phi[tid].y;
  }
}

__global__ void Normalization(cufftDoubleComplex *d_dfdc,
		              cufftDoubleComplex *d_dfdphi,
                         double one_by_nxnynz, int nx, int ny, int nz){

  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
  unsigned int tid = iz + nz * (iy + ix * ny);
     
  if(tid < nx * ny * nz){
     d_dfdc[tid].x *= one_by_nxnynz;
     d_dfdc[tid].y = 0.0;
     d_dfdphi[tid].x *= one_by_nxnynz;
     d_dfdphi[tid].y = 0.0;
  }
}

__global__ void Bounds(cufftDoubleComplex *d_comp, 
		       cufftDoubleComplex *d_dfdphi,
		       int *d_boundsFlag, int nx, int ny, int nz){

  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
  unsigned int tid = iz + nz * (iy + ix * ny);

  if(tid < nx * ny * nz){
    if(d_comp[tid].x < -0.4 || d_comp[tid].x > 1.4){
      *d_boundsFlag = 0;
      if (tid == 0) printf("Comp out of bounds\n");
    }
    /*if(d_dfdphi[tid].x < -0.4 || d_dfdphi[tid].x > 1.4){
      *d_boundsFlag = 0;
      if (tid == 0) printf("Phi out of bounds\n");
    }*/
  }
}

__global__ void MaxError(double *d_tempreal, cufftDoubleComplex *d_comp, 
			int nx, int ny, int nz){
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
  unsigned int tid = iz + nz * (iy + ix * ny);

  if(tid < nx * ny * nz)
    d_tempreal[tid] = fabs(d_tempreal[tid] - d_comp[tid].x);

}

__global__ void UpdateTempreal(double *d_tempreal, cufftDoubleComplex *d_comp, 
			int nx, int ny, int nz){
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
  unsigned int iz = threadIdx.z + blockIdx.z * blockDim.z;
  unsigned int tid = iz + nz * (iy + ix * ny);

  if(tid < nx * ny * nz){
    d_tempreal[tid] = d_comp[tid].x;
  }
}

