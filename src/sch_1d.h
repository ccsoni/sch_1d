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
  float tnow,  dtime, tend;

  float rho;  // dt/(dx)^2
  float hbar;

  float xmax, xmin;
  float delta_x;

  float mass;

  char model_name[MODEL_NAME_LENGTH];
  
  int32_t output_indx, noutput;
  float *output_timing;
};

#define SQR(x)  ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))

// square root of 2*M_PI
#define ROOT_2PI (1.772453851)
// 4-th root of 2*M_PI
#define QUAD_ROOT_2PI (1.583233487)

void init_run(struct run_param *, int, char**);
void setup_IC_free(float complex *, float, float, float, struct run_param *);
void evolve_3pnt(float complex *, struct run_param*, float);
void output_data(float complex *, float *, float *velc, struct run_param *);
