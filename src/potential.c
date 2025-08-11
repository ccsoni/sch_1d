#include "sch_1d.h"

void calc_pot(double *pot, struct run_param *tr)
{
  for(int ix=0;ix<tr->nmesh_x;ix++) pot[ix] = 0.0;
}
