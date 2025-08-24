#include "sch_1d.h"

void calc_pot(double *pot, struct run_param *tr)
{
  double omega = 0.25;
  for(int ix=0;ix<tr->nmesh_x;ix++) {
    double x = tr->xmin + ((double)ix+0.5)*tr->delta_x;
    pot[ix] = 0.5*SQR(omega)*x*x;
  }
}
