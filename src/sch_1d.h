#pragma once
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define MODEL_NAME_LENGTH (256)

struct run_param {
  int32_t nstep;
  int32_t nmesh_x, nmesh_v;
  double tnow,  dtime, tend;

  double rho;  // dt/(dx)^2
  double hbar;

  // phase space resolution where sigma_x*sigma_v = hbar/2
  double sigma_x, sigma_v;

  double xmax, xmin;
  double delta_x;

  double vmax, vmin;
  double delta_v;

  double mass;

  char model_name[MODEL_NAME_LENGTH];

  int32_t output_indx, noutput;
  double *output_timing;
};

#define SQR(x)  ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))

// square root of 2*M_PI
#define ROOT_2PI (1.772453851)
// 4-th root of 2*M_PI
#define QUAD_ROOT_2PI (1.583233487)

void init_run(struct run_param *, int, char**);
double complex coherent_wavefunc(double, double, double, double, double);
void setup_IC_point(double complex *, double, double, double, struct run_param *);
void setup_IC_expand(double complex *, double, struct run_param *);
void evolve_3pnt_free(double complex *, struct run_param*, double);
void evolve_5pnt_free(double complex *, struct run_param*, double);
void evolve_3pnt(double complex *, double *, struct run_param*, double);
void evolve_5pnt(double complex *, double *, struct run_param*, double);
void evolve_7pnt(double complex *, double *, struct run_param*, double);
void output_data(double complex *, double *, double *velc, struct run_param *);
void output_df(double *, struct run_param *);
double complex analytic_psi(double, double, double, double, double, double);
void calc_prob(double complex *, double *, struct run_param *);
void calc_velc(double complex *, double *, struct run_param *);
void calc_pot(double *, struct run_param *);
void calc_df(double complex *, double *, struct run_param *);
