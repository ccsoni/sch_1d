#include "sch_1d.h"

void output_data(double complex *psi, double *dens, double *velc,
		 struct run_param *this_run)
{
  static char output_filename[MODEL_NAME_LENGTH+10];

  if(this_run->tnow > this_run->output_timing[this_run->output_indx]) {
    sprintf(output_filename, "%s_%d.dat",
	    this_run->model_name, this_run->output_indx);
    FILE *fp = fopen(output_filename, "w");

    for(int32_t ix=0;ix<this_run->nmesh_x;ix++) {
      double x = this_run->xmin + ((double)ix+0.5)*this_run->delta_x;

      double complex analytic = analytic_psi(x, this_run->tnow,
					    0.5, 2.0*M_PI, 0.1, this_run->hbar);
      
      fprintf(fp, "%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
	      x, dens[ix], velc[ix], creal(psi[ix]), cimag(psi[ix]),
	      creal(analytic), cimag(analytic), SQR(cabs(analytic)));
    }

    fclose(fp);
    this_run->output_indx++;
  }
}
