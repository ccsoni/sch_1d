#include "sch_1d.h"



void calc_dens(float complex *psi, float *dens, struct run_param *tr)
{
  tr->mass=0.0;
  for(int32_t im=0;im<tr->nmesh_x;im++) {
    dens[im] = psi[im]*conj(psi[im]);
    tr->mass += dens[im]*tr->delta_x;
  }
}

void calc_velc(float complex *psi, float *velc, struct run_param *tr)
{

  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    float complex dpsi;

    int32_t ixp=(ix+1) % tr->nmesh_x;
    int32_t ixm=(ix-1+tr->nmesh_x) % tr->nmesh_x;
    dpsi = (psi[ixp]-psi[ixm])/(2.0*tr->delta_x);

    velc[ix] = tr->hbar*cimagf(dpsi*conjf(psi[ix]));
  }
}




int main(int argc, char **argv)
{
  struct run_param this_run = {.nmesh_x = 0};
  float complex *psi;
  float *dens, *velc, *DF;

  init_run(&this_run, argc, argv);

  psi = (float complex *) malloc(sizeof(float complex)*this_run.nmesh_x);
  dens = (float *) malloc(sizeof(float)*this_run.nmesh_x);
  velc = (float *) malloc(sizeof(float)*this_run.nmesh_x);

  float x_bar = 0.5;
  float v_bar = 2.0*M_PI;
  float sigma_x = 0.05;

  setup_IC_free(psi, x_bar, v_bar, sigma_x, &this_run);

  calc_dens(psi, dens, &this_run);
  calc_velc(psi, velc, &this_run);

  printf("# dt = %12.4e\n", this_run.dtime);

  while(this_run.tnow < this_run.tend) {
    
    evolve_3pnt(psi, &this_run, this_run.dtime);
    calc_dens(psi, dens, &this_run);
    calc_velc(psi, velc, &this_run);

    this_run.tnow += this_run.dtime;

    if(this_run.tnow > this_run.output_timing[this_run.output_indx]) {
      output_data(psi, dens, velc, &this_run);
    }

    
    
  }
 
#if 0
  for(int32_t ix=0;ix<this_run.nmesh_x;ix++) {
    float x = this_run.xmin + ((float)ix+0.5)*this_run.delta_x;
    printf("%14.6e %14.6e %14.6e %14.6e\n",
	   x, dens[ix], crealf(psi[ix]), cimagf(psi[ix]));
  }
#endif

}
