#include "sch_1d.h"

void calc_dens(double complex *psi, double *dens, struct run_param *tr)
{
  tr->mass=0.0;
  for(int32_t im=0;im<tr->nmesh_x;im++) {
    dens[im] = psi[im]*conj(psi[im]);
    tr->mass += dens[im]*tr->delta_x;
  }
}

void calc_velc(double complex *psi, double *velc, struct run_param *tr)
{

  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    double x = tr->xmin + ((double)ix+0.5)*tr->delta_x;
    int32_t ixp1=(ix+1) % tr->nmesh_x;
    int32_t ixm1=(ix-1+tr->nmesh_x) % tr->nmesh_x;
    int32_t ixp2=(ix+2) % tr->nmesh_x;
    int32_t ixm2=(ix-2+tr->nmesh_x) % tr->nmesh_x;    

#ifdef __SECOND_ORDER__
    double complex d_psi = (psi[ixp1]-psi[ixm1])/(2.0*tr->delta_x);
#else 
    double complex d_psi = (8.0*(psi[ixp1]-psi[ixm1])-(psi[ixp2]-psi[ixm2]))/(12.0*tr->delta_x);
#endif
    double dens = psi[ix]*conj(psi[ix]);

    velc[ix] = tr->hbar*cimag(d_psi*conj(psi[ix]))/dens;
  }
}

