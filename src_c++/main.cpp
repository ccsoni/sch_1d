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

  setup_IC_free_particle(psi_1d, 0.5, 2.0*M_PI, 0.02, this_run);

  DF = static_cast<double*>(std::aligned_alloc(64, sizeof(double)*this_run.nmesh_x*this_run.nmesh_v));

  calc_dens(psi_1d, dens, this_run);
  calc_velc(psi_1d, velc, this_run);
  calc_DF(DF, psi_1d, this_run);

  while(this_run.tnow < this_run.tend) {

    if(this_run.nstep % 100 == 0) {
      std::cout << "# t = " << std::scientific << std::setprecision(6) << this_run.tnow << std::endl;
    }
    
    evolve_5pnt(psi_1d, this_run, this_run.dtime);

    if(this_run.tnow > this_run.next_output_timing()) {
      calc_dens(psi_1d, dens, this_run);
      calc_velc(psi_1d, velc, this_run);
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
