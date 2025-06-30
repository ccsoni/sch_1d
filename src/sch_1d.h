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

#define MODEL_NAME_LENGTH (256)

struct run_param {
  int32_t nstep;
  int32_t nmesh_x;
  double tnow,  dtime, tend;

  double rho;  // dt/(dx)^2
  double hbar;

  double xmax, xmin;
  double delta_x;

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
void setup_IC_free(double complex *, double, double, double, struct run_param *);
void evolve_3pnt(double complex *, struct run_param*, double);
void evolve_5pnt(double complex *, struct run_param*, double);
void output_data(double complex *, double *, double *velc, struct run_param *);
double complex analytic_psi(double, double, double, double, double, double);
void calc_dens(double complex *, double *, struct run_param *);
void calc_velc(double complex *, double *, struct run_param *);
