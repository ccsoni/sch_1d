#include "sch_1d.h"
#include "prototype.h"

void output_data(complexd *psi, double *dens, double *velc, run_param & tr)
{
  static char output_filename[MODEL_NAME_LENGTH+10];

  std::snprintf(output_filename, sizeof(output_filename),
		"%s_%03d.dat", tr.model_name.c_str(), tr.output_indx);  

  FILE *fp = fopen(output_filename, "w");

  for(int32_t ix=0;ix<tr.nmesh_x;ix++) {
    double x = tr.xmin + ((double)ix+0.5)*tr.delta_x;

    complexd analytic = analytic_psi(x, tr.tnow,
				     0.0, 2.0*M_PI, 0.05, tr.hbar);

    fprintf(fp, "# tnow = %14.6e\n", tr.tnow);
      
    fprintf(fp, "%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
	    x, dens[ix], velc[ix], psi[ix].real(), psi[ix].imag(),
	    analytic.real(), analytic.imag(), std::norm(analytic));
  }
  fflush(fp);

  fclose(fp);
}

void output_DF(double *DF, run_param & tr)
{
  static char output_filename[MODEL_NAME_LENGTH+10];

  std::snprintf(output_filename, sizeof(output_filename),
		"%s_DF_%03d.dat", tr.model_name.c_str(), tr.output_indx);

  FILE *fp_DF = fopen(output_filename, "w");

  for(int32_t ix=0;ix<tr.nmesh_x;ix++) {
    for(int32_t iv=0;iv<tr.nmesh_v;iv++) {
      double x = tr.xmin + (static_cast<double>(ix)+0.5)*tr.delta_x;
      double v = tr.vmin + (static_cast<double>(iv)+0.5)*tr.delta_v;
      fprintf(fp_DF, "%14.6e %14.6e %14.6e\n",
	      x, v, DF[iv + tr.nmesh_v*ix]);
    }
    fprintf(fp_DF, "\n");
  }
 fclose(fp_DF);  
}
