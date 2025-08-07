#include "sch_1d.h"

double complex coherent_wavefunc(double x, double q, double v, double sigma_x,
				 double hbar)
{
  double complex wf = cexp(_Complex_I*v*x/hbar)*exp(-0.25*SQR((x-q)/sigma_x));
  double C= QUAD_ROOT_2PI*sqrt(sigma_x);

  return wf/C;
}

double sigma_norm(double t, double sigma, double hbar)
{
  double sig = 1.0 + SQR(hbar*t)/(4.0*QUAD(sigma));

  return sig;
}

double complex analytic_psi(double x, double t,
			    double q, double v, double sigma_x, double hbar)
{
  double sig = sigma_norm(t, sigma_x, hbar);
  double dx = x-q;

  double exp_arg_real = -SQR(x-q-v*t)/(4.0*SQR(sigma_x)*sig);
  double exp_arg_imag = SQR(x-q)*hbar*t/(8.0*QUAD(sigma_x)) + v*(x-q)/hbar - 0.5*SQR(v)*t/hbar;

  double complex exp_arg = exp_arg_real + _Complex_I*exp_arg_imag;
  
  double complex wf = cexp(exp_arg)/(QUAD_ROOT_2PI*sqrt(sigma_x));

  return wf;
}


// initial condition for free particles
void setup_IC_free(double complex *psi,
		   double x_bar, double v_bar, double sigma_x,
		   struct run_param *tr)
{
  assert(tr->nmesh_x != 0);
  
  tr->xmin = 0.0;
  tr->xmax = 1.0;
  tr->delta_x = (tr->xmax - tr->xmin)/tr->nmesh_x;

  tr->dtime = tr->rho*SQR(tr->delta_x);

  int32_t nmirror = 0.5/sigma_x;

#if 0
  for(int32_t im=0;im<tr->nmesh_x;im++) {
    double x = tr->xmin + ((double)im+0.5)*tr->delta_x;
    psi[im] = coherent_wavefunc(x, x_bar, v_bar, sigma_x, tr->hbar);
  }
#else
  for(int32_t im=0;im<tr->nmesh_x;im++) {
    psi[im] = 0.0;
    for(int32_t m=-nmirror;m<nmirror;m++) {
      double x = tr->xmin + ((double)im+0.5)*tr->delta_x
	+ (double)m*(tr->xmax-tr->xmin);
      psi[im] += coherent_wavefunc(x, x_bar, v_bar, sigma_x, tr->hbar);
    }
  }
#endif

  // range of velocity in phase space based on the Nyquist wavelength
  tr->vmin = -M_PI*tr->hbar/tr->delta_x;
  tr->vmax =  M_PI*tr->hbar/tr->delta_x;

  // mesh spacing in the velocity space is set to half of the one
  // obtained with the unceartainty principle
  tr->delta_v = 0.25*tr->hbar/tr->delta_x;

  tr->nmesh_v = (tr->vmax-tr->vmin)/tr->delta_v;
  
}

