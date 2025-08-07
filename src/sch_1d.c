#include "sch_1d.h"

int main(int argc, char **argv)
{
  struct run_param this_run = {.nmesh_x = 0};
  double complex *psi;
  double *dens, *velc, *df;

  init_run(&this_run, argc, argv);

  psi = (double complex *) malloc(sizeof(double complex)*this_run.nmesh_x);
  dens = (double *) malloc(sizeof(double)*this_run.nmesh_x);
  velc = (double *) malloc(sizeof(double)*this_run.nmesh_x);

  double x_bar = 0.5;
  double v_bar = 2.0*M_PI;
  double sigma_x = 0.1;

  setup_IC_free(psi, x_bar, v_bar, sigma_x, &this_run);

  df = (double *) malloc(sizeof(double)*this_run.nmesh_x*this_run.nmesh_v);  

  calc_dens(psi, dens, &this_run);
  calc_velc(psi, velc, &this_run);
  calc_df(psi, df, &this_run);

  printf("# dt = %12.4e\n", this_run.dtime);

  while(this_run.tnow < this_run.tend) {

#ifdef __SECOND_ORDER__
    evolve_5pnt(psi, &this_run, this_run.dtime);
#else
    evolve_5pnt(psi, &this_run, this_run.dtime);
#endif
    calc_dens(psi, dens, &this_run);
    calc_velc(psi, velc, &this_run);
    calc_df(psi, df, &this_run);

    this_run.tnow += this_run.dtime;

    if(this_run.tnow > this_run.output_timing[this_run.output_indx]) {
      output_data(psi, dens, velc, &this_run);
      output_df(df, &this_run);
      this_run.output_indx++;
    }
    
  }
 
  free(psi);
  free(dens);
  free(velc);
}
