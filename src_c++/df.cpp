#include "sch_1d.h"
#include "prototype.h"

// compute the phase space density at (x_, v_)
double calc_df_at(double x_, double v_, complexd *psi, run_param & tr)
{
  double kernel_norm = sqrt(2.0*M_PI*tr.hbar);

  // In computing the DF from the wave fuction, the sigma_x in the
  // Husimi integral is set to mesh spacing of the wave function.

  complexd phi_H = 0.0;
#pragma omp parallel for reduction(+:phi_H)
  for(int32_t ix=0;ix<tr.nmesh_x;ix++) {
    double x_int = tr.xmin + (static_cast<double>(ix)+0.5)*tr.delta_x;
    complexd W_kernel = coherent_wavefunc(x_int, x_, v_, tr.sigma_x, tr.hbar);

    phi_H += conj(W_kernel)/kernel_norm*psi[ix]*tr.delta_x;
  }

  return norm(phi_H);
}

void calc_DF(double *DF, complexd *psi, run_param & tr)
{
#pragma omp parallel for schedule(auto) collapse(2)
  for(int32_t ix=0;ix<tr.nmesh_x;ix++) {
    for(int32_t iv=0;iv<tr.nmesh_v;iv++) {
      double x_ = tr.xmin + (static_cast<double>(ix)+0.5)*tr.delta_x;
      double v_ = tr.vmin + (static_cast<double>(iv)+0.5)*tr.delta_v;
      // compute the phase space density at (x_, v_)
      DF[iv + tr.nmesh_v*ix] = calc_df_at(x_, v_, psi, tr);
    }
  }
}
