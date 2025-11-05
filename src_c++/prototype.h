#pragma once
#include "sch_1d.h"

void setup_IC_coherent_particle(complexd*, double, double, double, run_param &);
void setup_IC_expand(complexd*, double, run_param &);
void setup_IC_constvel(complexd*, double, run_param &);
void calc_prob(complexd*, double*, run_param &);
void calc_dens(complexd*, double*, run_param &);
void calc_velc(complexd*, double*, run_param &);
void calc_DF(double*, complexd*, run_param &);
void calc_pot(double*, double*, run_param &);
void calc_energy(double*, double*, run_param &);
void evolve_5pnt(complexd*, double*, run_param &, double);
void output_data(complexd*, double*, double*, run_param &);
void output_DF(double*, run_param &);
complexd analytic_psi(double, double, double, double, double, double);
complexd coherent_wavefunc(double, double, double, double, double);
