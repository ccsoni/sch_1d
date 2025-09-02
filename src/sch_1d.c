#include "sch_1d.h"

int main(int argc, char **argv)
{
  struct run_param this_run = {.nmesh_x = 0};
  double complex *psi;
  double *dens, *velc, *df;
  double *pot;

  init_run(&this_run, argc, argv);

  psi = (double complex *) malloc(sizeof(double complex)*this_run.nmesh_x);
  dens = (double *) malloc(sizeof(double)*this_run.nmesh_x);
  velc = (double *) malloc(sizeof(double)*this_run.nmesh_x);
  pot = (double *)malloc(sizeof(double)*this_run.nmesh_x);

  double x_bar = -0.5;
  double v_bar = 1.0;
  double sigma_x = 0.05;

  setup_IC_point(psi, x_bar, v_bar, sigma_x, &this_run);
  //setup_IC_expand(psi, v_bar, &this_run);

  df = (double *) malloc(sizeof(double)*this_run.nmesh_x*this_run.nmesh_v);  

  calc_df(psi, df, &this_run);
  calc_dens(psi, dens, &this_run);
  calc_velc(psi, velc, &this_run);

  calc_pot(pot, &this_run);

  printf("# dt = %12.4e\n", this_run.dtime);

  while(this_run.tnow < this_run.tend) {

    if(this_run.nstep % 100 == 0) printf("# nstep = %d / tnow = %14.6e / mass = %14.6e\n", this_run.nstep, this_run.tnow, this_run.mass);

#ifdef __3_PNT_APPROX__
    evolve_3pnt(psi, pot, &this_run, this_run.dtime);
#elif defined (__5_PNT_APPROX__)
    evolve_5pnt(psi, pot, &this_run, this_run.dtime);
#elif defined (__7_PNT_APPROX__)
    evolve_7pnt(psi, pot, &this_run, this_run.dtime);
#endif
    calc_dens(psi, dens, &this_run);
    calc_velc(psi, velc, &this_run);

    if(this_run.tnow > this_run.output_timing[this_run.output_indx]) {
      calc_df(psi, df, &this_run);
      output_data(psi, dens, velc, &this_run);
      output_df(df, &this_run);
      this_run.output_indx++;
    }

    this_run.tnow += this_run.dtime;
    this_run.nstep++;
    
  }
 
  free(psi);
  free(dens);
  free(velc);
}
