#ifndef _ODE_SOLVER_H_
#define _ODE_SOLVER_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include "../cell_models/noble_1962.h"
#include "../config/user_config.h"

struct ode_solver
{
    char model_name[MAX_FILENAME_SIZE];

    uint32_t num_ode_equations;
    uint32_t n_active_cells;

    double *sv;

    get_cell_model_data_fn *get_cell_model_data;
    set_ode_initial_conditions_cpu_fn *set_ode_initial_conditions_cpu;
    solve_model_ode_cpu_fn *solve_model_ode_cpu;

};

struct ode_solver* new_ode_solver ();
void configure_ode_solver (struct ode_solver *the_ode_solver, const uint32_t num_volumes);
void set_ode_initial_condition_for_all_volumes (struct ode_solver *the_ode_solver);
void print_state_vector (struct ode_solver *the_ode_solver);

#endif