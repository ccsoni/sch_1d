#include "sch_1d.h"
#include <fftw3.h>

void calc_prob(complexd *psi, double *prob, run_param & tr)
{
  tr.mass = 0.0;
#pragma omp parallel for schedule(auto) reduction(+:tr.mass)
  for(int32_t im=0;im<tr.nmesh_x;im++) {
    prob[im] = std::norm(psi[im]);
    tr.mass += prob[im]*tr.delta_x;
  }
}

void calc_dens(complexd *psi, double *dens, run_param & tr)
{
  static bool first_call = true;
  static fftw_plan forward, backward;
  static fftw_complex *dens_hat;

  if(first_call) {
    dens_hat =
      static_cast<fftw_complex*>(std::aligned_alloc(64, sizeof(fftw_complex)*(tr.nmesh_x/2 + 1)));

    forward = fftw_plan_dft_r2c_1d(tr.nmesh_x, dens, dens_hat, FFTW_ESTIMATE);
    backward = fftw_plan_dft_c2r_1d(tr.nmesh_x, dens_hat, dens, FFTW_ESTIMATE);

    first_call = false;
  }

  calc_prob(psi, dens, tr);

  fftw_execute(forward);

  // Gaussian filter in k-space
#pragma omp parallel for schedule(auto)
  for(int32_t ik=0;ik<tr.nmesh_x/2;ik++) {
    double kx = 2.0*M_PI*static_cast<double>(ik)/static_cast<double>(tr.nmesh_x);
    double wk = std::exp(-0.5*SQR(tr.sigma_x/tr.delta_x*kx));
    dens_hat[ik][0] *= wk;
    dens_hat[ik][1] *= wk;
  }

  fftw_execute(backward);

#pragma omp parallel for schedule(auto)
  for(int32_t ix=0;ix<tr.nmesh_x;ix++) dens[ix] /= tr.nmesh_x;

}


void calc_velc(complexd *psi, double *velc, run_param  & tr)
{

#pragma omp parallel for schedule(auto)
  for(int32_t ix=0;ix<tr.nmesh_x;ix++) {
    double x = tr.xmin + (static_cast<double>(ix)+0.5)*tr.delta_x;
    int32_t ixp1=(ix+1) % tr.nmesh_x;
    int32_t ixm1=(ix-1+tr.nmesh_x) % tr.nmesh_x;
    int32_t ixp2=(ix+2) % tr.nmesh_x;
    int32_t ixm2=(ix-2+tr.nmesh_x) % tr.nmesh_x;

#ifdef __SECOND_ORDER__
    complexd d_psi = (psi[ixp1]-psi[ixm1])/(2.0*tr.delta_x);
#else
    complexd d_psi = (8.0*(psi[ixp1]-psi[ixm1])-(psi[ixp2]-psi[ixm2]))/(12.0*tr.delta_x);
#endif
    double dens = norm(psi[ix]);

    velc[ix] = tr.hbar*imag(d_psi*conj(psi[ix]))/dens;
  }
}
