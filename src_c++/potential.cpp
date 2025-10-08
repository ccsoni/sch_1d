#include "sch_1d.h"

void calc_energy(double *DF, double *pot, run_param & tr)
{
  double K, W;
  K=0.0;
  W=0.0;

#pragma omp parallel for schedule(auto) reduction(+:K,W)
  for(int32_t ix=0;ix<tr.nmesh_x;ix++) {
    double v2_mean = 0.0;
    double dens = 0.0;
    for(int32_t iv=0;iv<tr.nmesh_v;iv++) {
      double v = tr.vmin + ((double)iv+0.5)*tr.delta_v;
      int32_t im = iv + tr.nmesh_v*ix;
      v2_mean += SQR(v)*DF[im]*tr.delta_v;
      dens += DF[im]*tr.delta_v;
    }
    v2_mean = v2_mean/dens;
    K += 0.5*v2_mean*dens*tr.delta_x;
    W += dens*pot[ix]*tr.delta_x;
  }

  tr.Kene = K;
  tr.Wene = W;
}

void calc_pot(double *pot, double *dens, run_param & tr)
{
  double omega = 2.0;;
  for(int32_t ix=0;ix<tr.nmesh_x;ix++) {
    double x = tr.xmin + (static_cast<double>(ix)+0.5)*tr.delta_x;
    pot[ix] = 0.5*SQR(omega*x);
    //    pot[ix] = 0.0;
  }
}
