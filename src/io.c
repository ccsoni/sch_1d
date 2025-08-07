#include "sch_1d.h"

void output_data(double complex *psi, double *dens, double *velc,
		 struct run_param *tr)
{
  static char output_filename[MODEL_NAME_LENGTH+10];

  sprintf(output_filename, "%s_%02d.dat",
	  tr->model_name, tr->output_indx);
  FILE *fp = fopen(output_filename, "w");

  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    double x = tr->xmin + ((double)ix+0.5)*tr->delta_x;

    double complex analytic = analytic_psi(x, tr->tnow,
					   0.5, 2.0*M_PI, 0.1, tr->hbar);
      
    fprintf(fp, "%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
	    x, dens[ix], velc[ix], creal(psi[ix]), cimag(psi[ix]),
	    creal(analytic), cimag(analytic), SQR(cabs(analytic)));
  }

  fclose(fp);
}
