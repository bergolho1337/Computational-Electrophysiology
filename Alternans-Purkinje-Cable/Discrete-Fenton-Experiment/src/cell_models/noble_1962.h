#ifndef _NOBLE_1962_H_
#define _NOBLE_1962_H_

#include <iostream>
#include "model_common.h"

#define NEQ 4
#define INITIAL_V (-75.5344986658)

static const double stim_current = 1000.0;
static const double stim_start = 0.0;
static const double stim_duration = 2.0;
static const double start_period = 300.0;
static const double end_period = 300.0;
static const double period_step = 100.0;
static const int n_cycles = 10;
static const int id_limit = 20;

void init_cell_model_data (uint32_t *num_equations);
void set_model_initial_conditions_cpu (double *sv);
void solve_model_odes_cpu (const double dt, double *sv, double *vstar,\
                            double *stim_currents, const uint32_t num_cells_to_solve);
void solve_model_ode_cpu(double dt, double *sv, double vstar, double stim_current);
void RHS_cpu(const double *sv, const double vstar, double *rDY_, double stim_current);

#endif