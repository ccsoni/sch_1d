#include "sch_1d.h"

int main(int argc, char **argv)
{
  struct run_param this_run = {.nmesh_x = 0};
  double complex *psi;
  double *dens, *prob, *velc, *df;
  double *pot;

  init_run(&this_run, argc, argv);

  psi = (double complex *) malloc(sizeof(double complex)*this_run.nmesh_x);
  dens = (double *) malloc(sizeof(double)*this_run.nmesh_x);
  prob = (double *) malloc(sizeof(double)*this_run.nmesh_x);
  velc = (double *) malloc(sizeof(double)*this_run.nmesh_x);
  pot = (double *)malloc(sizeof(double)*this_run.nmesh_x);

  double x_bar = 0.0;
  double v_bar = 1.0;
  double sigma_x = 0.05;
  setup_IC_coherent_particle(psi, x_bar, v_bar, sigma_x, &this_run);

  calc_dens(psi, dens, &this_run);

  calc_pot(pot, dens, &this_run);

  for(int32_t ix=1;ix<this_run.nmesh_x-1;ix++) {
    double x = this_run.xmin + ((double)ix+0.5)*this_run.delta_x;
    int32_t im = (ix - 1 + this_run.nmesh_x)%this_run.nmesh_x;
    int32_t ip = (ix + 1)%this_run.nmesh_x;
    printf("%14.6e %14.6e %14.6e %14.6e\n", x, dens[ix], pot[ix],
	   (pot[ip]-2.0*pot[ix]+pot[im])/(4.0*M_PI*SQR(this_run.delta_x)));
  }
}
