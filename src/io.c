#include "sch_1d.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>

#define MAXPATHLEN (1024)

void make_directory(char *directory_name)
{
  mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;

  static char cwd_path[MAXPATHLEN];
  char *path;

  path = getcwd(cwd_path, sizeof(cwd_path));

  strcat(cwd_path,"/");
  strcat(cwd_path, directory_name);

  mkdir(cwd_path, mode);
}

void output_data(double complex *psi, double *prob, double *velc,
		 struct run_param *tr)
{
  static char output_filename[2*MODEL_NAME_LENGTH+10];

  make_directory(tr->model_name);

  sprintf(output_filename, "%s/%s_%03d.dat", tr->model_name,
	  tr->model_name, tr->output_indx);
  FILE *fp = fopen(output_filename, "w");

  fprintf(fp, "%d\n", tr->nmesh_x);
  fprintf(fp, "%14.6e\n", tr->hbar);
  fprintf(fp, "%14.6e\n", tr->tnow);
  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    double x = tr->xmin + ((double)ix+0.5)*tr->delta_x;

    double complex analytic = analytic_psi(x, tr->tnow,
					   -0.5, 1.0, 0.05, tr->hbar);

    fprintf(fp, "%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
	    x, prob[ix], velc[ix], creal(psi[ix]), cimag(psi[ix]),
	    creal(analytic), cimag(analytic), SQR(cabs(analytic)));
  }

  fclose(fp);
}

void output_df(double *df, struct run_param *tr)
{
  static char output_filename[2*MODEL_NAME_LENGTH+10];
  FILE *fp_DF;

  make_directory(tr->model_name);

  sprintf(output_filename, "%s/%s_DF_%03d.dat",
	  tr->model_name,
          tr->model_name, tr->output_indx);

  fp_DF = fopen(output_filename, "w");

  fprintf(fp_DF, "%d %d\n", tr->nmesh_x, tr->nmesh_v);
  fprintf(fp_DF, "%14.6e\n", tr->hbar);
  fprintf(fp_DF, "%14.6e\n", tr->tnow);

  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    for(int32_t iv=0;iv<tr->nmesh_v;iv++) {
      double x = tr->xmin + ((double)ix+0.5)*tr->delta_x;
      double v = tr->vmin + ((double)iv+0.5)*tr->delta_v;
      fprintf(fp_DF, "%14.6e %14.6e %14.6e\n", x, v, df[iv + tr->nmesh_v*ix]);
    }
    fprintf(fp_DF, "\n");
  }
  fclose(fp_DF);
}
