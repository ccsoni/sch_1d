#include "sch_1d.h"

complexd coherent_wavefunc(double x, double q, double v, double sigma_x,
			   double hbar)
{
  complexd wf = std::polar(1.0, v*(x-q)/hbar)*exp(-0.25*SQR((x-q)/sigma_x));
  double C= QUAD_ROOT_2PI*sqrt(sigma_x);

  return wf/C;
}

double sigma_norm(double t, double sigma, double hbar)
{
  double sig = 1.0 + SQR(hbar*t)/(4.0*QUAD(sigma));

  return sig;
}

complexd analytic_psi(double x, double t,
		      double q, double v, double sigma_x, double hbar)
{
  double sig = sigma_norm(t, sigma_x, hbar);
  double dx = x-q;

  double exp_arg_real = -SQR(x-q-v*t)/(4.0*SQR(sigma_x)*sig);
  double exp_arg_imag = (SQR(x-q)*hbar*t/(8.0*QUAD(sigma_x)) + v*(x-q)/hbar - 0.5*SQR(v)*t/hbar)/sig;

  complexd exp_arg = complexd(exp_arg_real, exp_arg_imag);
  complexd fact = complexd(1.0, hbar*t/(2.0*SQR(sigma_x)));

  complexd wf = exp(exp_arg)/(QUAD_ROOT_2PI*sqrt(sigma_x))/sqrt(fact);

  return wf;
}

void setup_IC_free_particle(complexd *psi, double x_, double v_,
			    double sigma_x, run_param & tr)
{
  tr.tnow = 0.0;
  tr.nstep = 0;

  tr.xmin =  -1.0;
  tr.xmax =  1.0;

  tr.delta_x = (tr.xmax - tr.xmin)/static_cast<double>(tr.nmesh_x);
  tr.dtime = tr.rho*SQR(tr.delta_x);

  for(int32_t im=0;im<tr.nmesh_x;im++) {
    double x = tr.xmin + (static_cast<double>(im)+0.5)*tr.delta_x;
    psi[im] = coherent_wavefunc(x, x_, v_, sigma_x, tr.hbar);
  }

  // velocity range in the phase space based on the Nyquist wavelength
  tr.vmin = -M_PI*tr.hbar/tr.delta_x;
  tr.vmax =  M_PI*tr.hbar/tr.delta_x;

  // spatial resolution in the phase space
  tr.sigma_x = 4.0*tr.delta_x;
  tr.sigma_v = 0.5*tr.hbar/tr.sigma_x;

  // mesh spacing in the velocity space is set ot 1/4 of the one
  // obtained with the unceartainty principle
  tr.delta_v = tr.sigma_v/4.0;
  tr.nmesh_v = (tr.vmax-tr.vmin)/tr.delta_v;

  assert(fabs(v_) < tr.vmax);
}

// IC with v = A*x
void setup_IC_expand(complexd *psi, double expand_coeff, run_param & tr)
{
  assert(tr.nmesh_x != 0);

  tr.nstep = 0;

  tr.xmin = -1.0;
  tr.xmax =  1.0;
  tr.delta_x = (tr.xmax - tr.xmin)/tr.nmesh_x;

  tr.dtime = tr.rho*SQR(tr.delta_x);

  // density profile
  // rho(x) = rho_0*(1-|x|/xcut)
  tr.mass = 1.0;
  double xcut = 0.2;  // cut off distance beyond which the density is zero.
  double dens_0 = 1.0/(2.0*xcut);

  for(int32_t ix=0;ix<tr.nmesh_x;ix++) {
    double x = tr.xmin + (static_cast<double>(ix)+0.5)*tr.delta_x;
    double dens;
    if(fabs(x) < xcut) {
      dens = dens_0*(1.0-fabs(x)/xcut);
    }else{
      dens = 0.0;
    }

    double theta = 0.5*expand_coeff*SQR(x);
    complexd arg = complexd(0.0, theta/tr.hbar);
    psi[ix] = sqrt(dens)*exp(arg);
  }

  // velocity range in the phase space based on the Nyquist wavelength
  tr.vmin = -M_PI*tr.hbar/tr.delta_x;
  tr.vmax =  M_PI*tr.hbar/tr.delta_x;

  // spatial resolution in the phase space
  tr.sigma_x = 4.0*tr.delta_x;
  tr.sigma_v = 0.5*tr.hbar/tr.sigma_x;

  // mesh spacing in the velocity space is set ot 1/4 of the one
  // obtained with the unceartainty principle
  tr.delta_v = tr.sigma_v/4.0;
  tr.nmesh_v = (tr.vmax-tr.vmin)/tr.delta_v;

  assert(fabs(expand_coeff*xcut) < tr.vmax);
}
