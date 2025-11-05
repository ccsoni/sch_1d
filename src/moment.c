#include "sch_1d.h"
#include "fftw3.h"

void calc_prob(double complex *psi, double *prob, struct run_param *tr)
{
  tr->mass=0.0;
  for(int32_t im=0;im<tr->nmesh_x;im++) {
    prob[im] = psi[im]*conj(psi[im]);
    tr->mass += prob[im]*tr->delta_x;
  }
}

void calc_dens(double complex *psi, double *dens, struct run_param *tr)
{
  static bool first_call = true;
  static fftw_plan forward, backward;
  static fftw_complex *dens_hat;

  if(first_call) {
    dens_hat =
      (fftw_complex *)aligned_alloc(64, sizeof(fftw_complex)*(tr->nmesh_x/2 + 1));

    forward = fftw_plan_dft_r2c_1d(tr->nmesh_x, dens, dens_hat, FFTW_ESTIMATE);
    backward = fftw_plan_dft_c2r_1d(tr->nmesh_x, dens_hat, dens, FFTW_ESTIMATE);

    first_call = false;
  }

  calc_prob(psi, dens, tr);

  return;

  fftw_execute(forward);

#if 0
  // Gaussian filter in k-space
#pragma omp parallel for schedule(auto)
  for(int32_t ik=0;ik<tr->nmesh_x/2;ik++) {
    double kx = 2.0*M_PI*(double)ik/(double)tr->nmesh_x;
    double wk = exp(-0.5*SQR(tr->sigma_x/tr->delta_x*kx));
    dens_hat[ik] *= wk;
  }
#endif

  fftw_execute(backward);

#pragma omp parallel for schedule(auto)
  for(int32_t ix=0;ix<tr->nmesh_x;ix++) dens[ix] /= tr->nmesh_x;
}

void calc_velc(double complex *psi, double *velc, struct run_param *tr)
{

  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    double x = tr->xmin + ((double)ix+0.5)*tr->delta_x;
    int32_t ixp1=(ix+1) % tr->nmesh_x;
    int32_t ixm1=(ix-1+tr->nmesh_x) % tr->nmesh_x;
    int32_t ixp2=(ix+2) % tr->nmesh_x;
    int32_t ixm2=(ix-2+tr->nmesh_x) % tr->nmesh_x;

#ifdef __SECOND_ORDER__
    double complex d_psi = (psi[ixp1]-psi[ixm1])/(2.0*tr->delta_x);
#else
    double complex d_psi = (8.0*(psi[ixp1]-psi[ixm1])-(psi[ixp2]-psi[ixm2]))/(12.0*tr->delta_x);
#endif
    double prob = psi[ix]*conj(psi[ix]);

    velc[ix] = tr->hbar*cimag(d_psi*conj(psi[ix]))/prob;
  }
}
