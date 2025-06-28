#include "sch_1d.h"

void output_data(float complex *psi, float *dens, float *velc,
		 struct run_param *this_run)
{
  static char output_filename[MODEL_NAME_LENGTH+10];

  if(this_run->tnow > this_run->output_timing[this_run->output_indx]) {
    sprintf(output_filename, "%s_%d.dat",
	    this_run->model_name, this_run->output_indx);
    FILE *fp = fopen(output_filename, "w");

    for(int32_t ix=0;ix<this_run->nmesh_x;ix++) {
      float x = this_run->xmin + ((float)ix+0.5)*this_run->delta_x;
      fprintf(fp, "%14.6e %14.6e %14.6e %14.6e %14.6e\n",
	      x, dens[ix], velc[ix], crealf(psi[ix]), cimagf(psi[ix]));
    }

    fclose(fp);
    this_run->output_indx++;
  }
}
