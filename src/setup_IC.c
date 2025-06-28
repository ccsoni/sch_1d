#include "sch_1d.h"

float complex coherent_wavefunc(float x, float q, float v, float sigma_x,
				float hbar)
{
  float complex wf = cexp(-_Complex_I*v*x/hbar)*exp(-0.25*SQR((x-q)/sigma_x));
  float C= QUAD_ROOT_2PI*sqrt(sigma_x);

  return wf/C;
}


// initial condition for free particles
void setup_IC_free(float complex *psi,
		   float x_bar, float v_bar, float sigma_x,
		   struct run_param *tr)
{
  assert(tr->nmesh_x != 0);
  
  tr->xmin = 0.0;
  tr->xmax = 1.0;
  tr->delta_x = (tr->xmax - tr->xmin)/tr->nmesh_x;

  tr->dtime = tr->rho*SQR(tr->delta_x);

  int32_t nmirror = 0.5/sigma_x;

  for(int32_t im=0;im<tr->nmesh_x;im++) {
    psi[im] = 0.0;
    for(int32_t m=-nmirror;m<nmirror;m++) {
      float x = tr->xmin + ((float)im+0.5)*tr->delta_x
	+ (float)m*(tr->xmax-tr->xmin);
      psi[im] += coherent_wavefunc(x, x_bar, v_bar, sigma_x, tr->hbar);
    }
  }
}
