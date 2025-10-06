#include "sch_1d.h"

// compute the phase space density at (x_, v_)
double calc_df_at(double x_, double v_, double complex *psi, struct run_param *tr)
{
  double kernel_norm = sqrt(2.0*M_PI*tr->hbar);

  // In computing the DF from the wave fuction, the sigma_x in the
  // Husimi integral is set to mesh spacing of the wave function.

  complex double phi_H = 0.0;
  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    double x_int = tr->xmin + ((double)ix+0.5)*tr->delta_x;
    double complex W_kernel = coherent_wavefunc(x_int, x_, v_, 4.0*tr->delta_x, tr->hbar);

    phi_H += conj(W_kernel)/kernel_norm*psi[ix]*tr->delta_x;
  }

  return phi_H*conj(phi_H);
}

void calc_df(double complex *psi, double *DF, double *dens,
	     struct run_param *tr)
{
#pragma omp parallel for schedule(auto)
  for(int32_t ix=0;ix<tr->nmesh_x;ix++) dens[ix] = 0.0;

#pragma omp parallel for schedule(auto) collapse(2)
  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    for(int32_t iv=0;iv<tr->nmesh_v;iv++) {
      // compute the phase space density at (x_, v_)
      double x_ = tr->xmin + ((double)ix+0.5)*tr->delta_x;
      double v_ = tr->vmin + ((double)iv+0.5)*tr->delta_v;
      DF[iv + tr->nmesh_v*ix] = calc_df_at(x_, v_, psi, tr);
      dens[ix] += DF[iv + tr->nmesh_v*ix]*tr->delta_v;
    }
  }
}
