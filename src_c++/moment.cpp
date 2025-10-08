#include "sch_1d.h"

void calc_prob(complexd *psi, double *prob, run_param & tr)
{
  tr.mass = 0.0;
  for(int32_t im=0;im<tr.nmesh_x;im++) {
    prob[im] = std::norm(psi[im]);
    tr.mass += prob[im]*tr.delta_x;
  }
}


void calc_velc(complexd *psi, double *velc, run_param  & tr)
{

  for(int32_t ix=0;ix<tr.nmesh_x;ix++) {
    double x = tr.xmin + (static_cast<double>(ix)+0.5)*tr.delta_x;
    int32_t ixp1=(ix+1) % tr.nmesh_x;
    int32_t ixm1=(ix-1+tr.nmesh_x) % tr.nmesh_x;
    int32_t ixp2=(ix+2) % tr.nmesh_x;
    int32_t ixm2=(ix-2+tr.nmesh_x) % tr.nmesh_x;

#ifdef __SECOND_ORDER__
    complexd d_psi = (psi[ixp1]-psi[ixm1])/(2.0*tr.delta_x);
#else
    complexd d_psi = (8.0*(psi[ixp1]-psi[ixm1])-(psi[ixp2]-psi[ixm2]))/(12.0*tr.delta_x);
#endif
    double dens = norm(psi[ix]);

    velc[ix] = tr.hbar*imag(d_psi*conj(psi[ix]))/dens;
  }
}
