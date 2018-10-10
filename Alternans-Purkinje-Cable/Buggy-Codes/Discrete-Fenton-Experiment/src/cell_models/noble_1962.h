#ifndef _NOBLE_1962_H_
#define _NOBLE_1962_H_

#include <iostream>
#include "model_common.h"
#include "../cell/cell.h"

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

//___________________________________________________________________________
//Parameters (miliseconds)
static const double Cm = 12.0;                                 // (microF)
static const double g_na_max = 400.0;                       // (microS)
static const double E_na = 40.0;                               // (millivolt)
static const double g_L = 0.075;                                // (microS)
static const double E_L = -60.0;                               // (millivolt)

void init_cell_model_data (uint32_t *num_equations);
void set_model_initial_conditions_cpu (double *sv);
void solve_model_odes_cpu (struct cell_data *volumes, double *stim_currents,\
                           const double dt ,const uint32_t num_cells_to_solve);
double dvdt (double *y_star, double stim_current);
double dmdt (double *y_old);
double dhdt (double *y_old);
double dndt (double *y_old);

#endif