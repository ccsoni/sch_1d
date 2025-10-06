#include "sch_1d.h"

void init_run(struct run_param *tr, int argc, char **argv)
{
  // default setting
  double dt_output = 0.1;
  tr->nmesh_x = 128;
  tr->rho = 1.0;
  tr->hbar = 1.0e-1;

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
      tr->nmesh_x = atof(optarg);
      fprintf(stderr, "# nmesh_x = %d\n", tr->nmesh_x);
      break;
    case 'r':
      tr->rho = atof(optarg);
      fprintf(stderr, "# dt/(dx)^2 = %12.4e\n", tr->rho);
      break;
    case 'T':
      tr->tend = atof(optarg);
      Tflag = true;
      fprintf(stderr, "# t_end = %12.4e\n", tr->tend);
      break;
    case 't':
      dt_output = atof(optarg);
      tflag = true;
      fprintf(stderr, "# dt_output = %12.4e\n", dt_output);
      break;
    case 'm':
      sprintf(tr->model_name, "%s", optarg);
      mflag = true;
      fprintf(stderr, "# model_name = %s\n", tr->model_name);
      break;
    case 'd':
      tr->hbar = atof(optarg);
      fprintf(stderr, "# hbar = %12.4e\n", tr->hbar);
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

  tr->output_indx = 0;
  tr->tnow = 0.0f;
  tr->noutput = (tr->tend-tr->tnow)/dt_output+1;
  tr->output_timing = (double *)malloc(sizeof(double)*tr->noutput);

  fprintf(stderr, "# nmesh_x : %d\n", tr->nmesh_x);
  fprintf(stderr, "# dt/dx^2 : %12.4e\n", tr->rho);
  fprintf(stderr, "# hbar    : %12.4e\n", tr->hbar);

  for(int i=0;i<tr->noutput-1;i++) {
    tr->output_timing[i] = tr->tnow + (i+1)*dt_output;
    fprintf(stderr, "%d %12.4e \n", i, tr->output_timing[i]);
  }

  tr->output_timing[tr->noutput-1] = tr->tend;
}
