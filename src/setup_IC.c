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

  double complex exp_arg = exp_arg_real + _Complex_I*exp_arg_imag/sig;
  double complex fact = 1.0/csqrt(1.0+_Complex_I*hbar*t/(2.0*SQR(sigma_x)));
  
  double complex wf = cexp(exp_arg)/(QUAD_ROOT_2PI*sqrt(sigma_x))*fact;

  return wf;
}


// initial condition for point mass like particles
void setup_IC_point(double complex *psi,
		   double x_bar, double v_bar, double sigma_x,
		   struct run_param *tr)
{
  assert(tr->nmesh_x != 0);
  
  tr->xmin = -1.0;
  tr->xmax = 1.0;
  tr->delta_x = (tr->xmax - tr->xmin)/tr->nmesh_x;

  tr->dtime = tr->rho*SQR(tr->delta_x);

  int32_t nmirror = 0.5/sigma_x;

#if 1
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

  // spatial resolution in the phase space
  tr->sigma_x = 4.0*tr->delta_x;
  tr->sigma_v = 0.5*tr->hbar/tr->sigma_x;
  
  // mesh spacing in the velocity space is set to half of the one
  // obtained with the unceartainty principle
  tr->delta_v = tr->sigma_v/4.0;

  tr->nmesh_v = (tr->vmax-tr->vmin)/tr->delta_v;
  
}

// IC with v = A*x and \rho(x) = 
void setup_IC_expand(double complex *psi, double expand_coeff, struct run_param *tr)
{
  assert(tr->nmesh_x != 0);
  
  tr->xmin = -1.0;
  tr->xmax = 1.0;
  tr->delta_x = (tr->xmax - tr->xmin)/tr->nmesh_x;

  tr->dtime = tr->rho*SQR(tr->delta_x);
  
  //density profile
  tr->mass = 1.0;
  double sigma_x = 0.3;
  double x0 = 0.0;
  double xcut=0.5;

  for(int32_t ix=0;ix<tr->nmesh_x;ix++) {
    double x = tr->xmin + ((double)ix+0.5)*tr->delta_x;

    double dens = (fabs(x) < xcut) ? 1.0-fabs(x)/xcut : 0.0;
 
    double theta = 0.5*expand_coeff*SQR(x);
    psi[ix] = sqrt(dens)*cexp(_Complex_I*theta/tr->hbar);
  }

  // range of velocity in phase space based on the Nyquist wavelength
  tr->vmin = -M_PI*tr->hbar/tr->delta_x;
  tr->vmax =  M_PI*tr->hbar/tr->delta_x;

  // spatial resolution in the phase space
  tr->sigma_x = 4.0*tr->delta_x;
  tr->sigma_v = 0.5*tr->hbar/tr->sigma_x;

  // mesh spacing in the velocity space is set to half of the one
  // obtained with the unceartainty principle
  tr->delta_v = tr->sigma_v/4.0;

  tr->nmesh_v = (tr->vmax-tr->vmin)/tr->delta_v;

}
