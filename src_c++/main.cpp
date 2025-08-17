#include "sch_1d.h"
#include "prototype.h"

int main(int argc, char **argv)
{
  run_param this_run;

  complexd *psi_1d;
  double *dens, *velc, *DF;

  this_run.init_run(argc, argv);

  psi_1d = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*this_run.nmesh_x));
  dens = static_cast<double*>(std::aligned_alloc(64, sizeof(double)*this_run.nmesh_x));
  velc = static_cast<double*>(std::aligned_alloc(64, sizeof(double)*this_run.nmesh_x));

  double x_ = 0.0;
  double v_ = 2.0*M_PI;
  double sigma_x_ = 0.05;

  setup_IC_free_particle(psi_1d, x_, v_, sigma_x_, this_run);

  DF = static_cast<double*>(std::aligned_alloc(64, sizeof(double)*this_run.nmesh_x*this_run.nmesh_v));

  calc_dens(psi_1d, dens, this_run);
  calc_velc(psi_1d, velc, this_run);
  calc_DF(DF, psi_1d, this_run);

  while(this_run.tnow < this_run.tend) {

    if(this_run.nstep % 100 == 0) printf("# nstep = %d / tnow = %14.6e / mass = %14.6e\n", this_run.nstep, this_run.tnow, this_run.mass);

#ifdef __SECOND_ORDER__
    evolve_3pnt(psi_1d, this_run, this_run.dtime);
#else
    evolve_5pnt(psi_1d, this_run, this_run.dtime);
#endif

    calc_dens(psi_1d, dens, this_run);
    calc_velc(psi_1d, velc, this_run);

    if(this_run.tnow > this_run.next_output_timing()) {
      calc_DF(DF, psi_1d, this_run);
      
      output_data(psi_1d, dens, velc, this_run);
      output_DF(DF, this_run);
      this_run.output_indx++;
    }

    this_run.tnow += this_run.dtime;
    this_run.nstep++;
  }
  
  std::free(psi_1d);
}
