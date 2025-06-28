#include "sch_1d.h"

void init_run(struct run_param *this_run, int argc, char **argv)
{
  // default setting
  float dt_output = 0.1;
  this_run->nmesh_x = 128;
  this_run->rho = 1.0;
  this_run->hbar = 1.0e-2;

  bool mflag, tflag, Tflag;

  mflag = tflag = Tflag = false;

  int32_t c = getopt(argc, argv, "N:r:T:t:m:d:h");

  do {
    switch(c) {
    case 'h':
      fprintf(stderr, "Mandatory options:\n -T t_end\n -t dt_output\n -m model_name\n");
      fprintf(stderr, "Additional options:\n -d h_bar\n -r dt/(dx^2)\n");
      exit(EXIT_FAILURE);
    case 'N':
      this_run->nmesh_x = atof(optarg);
      fprintf(stderr, "# nmehs_x = %d\n", this_run->nmesh_x);
      break;
    case 'r':
      this_run->rho = atof(optarg);
      fprintf(stderr, "# dt/(dx)^2 = %12.4e\n", this_run->rho);
      break;
    case 'T':
      this_run->tend = atof(optarg);
      Tflag = true;
      fprintf(stderr, "# t_end = %12.4e\n", this_run->tend);
      break;
    case 't':
      dt_output = atof(optarg);
      tflag = true;
      fprintf(stderr, "# dt_output = %12.4e\n", dt_output);
      break;
    case 'm':
      sprintf(this_run->model_name, "%s", optarg);
      mflag = true;
      fprintf(stderr, "# model_name = %s\n", this_run->model_name);
      break;
    case 'd':
      this_run->hbar = atof(optarg);
      fprintf(stderr, "# hbar = %12.4e\n", this_run->hbar);
      break;
    default:
      fprintf(stderr,
	      "Usage: %s -T <end_time> -t <output_interval> -m <model_name>\n",
	      argv[0]);
      exit(EXIT_FAILURE);
    }
  } while((c=getopt(argc, argv, "N:r:T:t:m:d:h")) != -1);

  if(((Tflag && tflag) && mflag) == false) {
    fprintf(stderr,
	    "Usage: %s -T <end_time> -t <output_interval> -m <model_name>\n",
	    argv[0]);
    exit(EXIT_FAILURE);
  }

  this_run->output_indx = 0;
  this_run->tnow = 0.0f;
  this_run->dtime = this_run->rho*SQR(this_run->delta_x);
  this_run->noutput = (this_run->tend-this_run->tnow)/dt_output+1;
  this_run->output_timing = (float *)malloc(sizeof(float)*this_run->noutput);

  fprintf(stderr, "# nmesh_x : %d\n", this_run->nmesh_x);
  fprintf(stderr, "# dt/dx^2 : %12.4e\n", this_run->rho);
  fprintf(stderr, "# hbar    : %12.4e\n", this_run->hbar);

  for(int i=0;i<this_run->noutput-1;i++) {
    this_run->output_timing[i] = this_run->tnow + (i+1)*dt_output;
    fprintf(stderr, "%d %12.4e \n", i, this_run->output_timing[i]);
  }

  this_run->output_timing[this_run->noutput-1] = this_run->tend;
}
