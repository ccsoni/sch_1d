#include "sch_1d.h"
#include "prototype.h"

#pragma omp declare reduction(+: complexd: omp_out += omp_in) \
                    initializer(omp_priv = complexd(0.0, 0.0))

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
    //complexd W_kernel = coherent_wavefunc(x_int, x_, v_, tr.sigma_x, tr.hbar);
    complexd W_kernel = coherent_wavefunc(x_int, x_, v_, tr.sigma_x_array[ix], tr.hbar);

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


void calc_density_from_df(double *density, double *DF, run_param & tr)
{
#pragma omp parallel for schedule(auto)
  for(int32_t ix=0; ix<tr.nmesh_x; ix++) {
    double sum_f = 0.0;
    for(int32_t iv=0; iv<tr.nmesh_v; iv++) {
      sum_f += DF[iv + tr.nmesh_v * ix];
    }
    density[ix] = sum_f * tr.delta_v;
  }
}


void calc_velocity_from_df(double *velocity, double *density, double *DF, run_param & tr)
{
  const double eps = 1e-20;
  for(int32_t ix=0; ix<tr.nmesh_x; ix++) {
    double sum_vf = 0.0;
    for(int32_t iv=0; iv<tr.nmesh_v; iv++) {
      double v_ = tr.vmin + (static_cast<double>(iv) + 0.5) * tr.delta_v;
      sum_vf += v_ * DF[iv + tr.nmesh_v * ix];
    }
    velocity[ix] = sum_vf * tr.delta_v;
    //double flux_j = sum_vf * tr.delta_v;
    //velocity[ix] = flux_j / (density[ix] + eps);
  }
}
