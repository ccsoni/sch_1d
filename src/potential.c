#include "sch_1d.h"

void calc_energy(double *DF, double *pot, struct run_param *tr)
{
  double K, W;
  K=0.0;
  W=0.0;

#pragma omp parallel for schedule(auto) reduction(+:K,W)
  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    double v2_mean = 0.0;
    double dens = 0.0;
    for(int32_t iv=0;iv<tr->nmesh_v;iv++) {
      double v = tr->vmin + ((double)iv+0.5)*tr->delta_v;
      int32_t im = iv + tr->nmesh_v*ix;
      v2_mean += SQR(v)*DF[im]*tr->delta_v;
      dens += DF[im]*tr->delta_v;
    }
    v2_mean = v2_mean/dens;
    K += 0.5*v2_mean*dens*tr->delta_x;
    W += dens*pot[ix]*tr->delta_x;
  }

  tr->Kene = K;
  tr->Wene = W;
}

void calc_pot(double *pot, double *dens, struct run_param *tr)
{
#ifdef __SELF_GRAV__
#pragma omp parallel for schedule(auto)
  for(int32_t im=0;im<tr->nmesh_x;im++) {
    pot[im] = 0.0;
    double xi = tr->xmin + ((double)im+0.5)*tr->delta_x;
    for(int32_t jm=0;jm<im;jm++) {
      double xj = tr->xmin + ((double)jm+0.5)*tr->delta_x;
      pot[im] += 2.0*M_PI*(xi-xj)*dens[jm]*tr->delta_x;
    }
    for(int32_t jm=im;jm<tr->nmesh_x;jm++) {
      double xj = tr->xmin + ((double)jm+0.5)*tr->delta_x;
      pot[im] -= 2.0*M_PI*(xi-xj)*dens[jm]*tr->delta_x;
    }
  }
#elif defined __GRAV__
  // gravitational potential of Gaussian mass distribution of mass M_center centered at x=0
  double M_center = 1.0;
  double sigma = tr->sigma_x;
  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    double x = tr->xmin + ((double)ix+0.5)*tr->delta_x;
    //    pot[ix] = 2.0*M_PI*M_center*fabs(x);
    pot[ix] = 2.0*M_PI*M_center*(x*erf(x/sqrt(2.0*SQR(sigma))) + sqrt(2.0/M_PI)*sigma*exp(-0.5*SQR(x/sigma)));
  }
#elif defined __HARMONIC__
  double omega = 2.0;
  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    double x = tr->xmin + ((double)ix+0.5)*tr->delta_x;
    pot[ix] = 0.5*SQR(omega)*x*x;
  }
#else
  for(int ix=0;ix<tr->nmesh_x;ix++) pot[ix] = 0.0;
#endif
}
