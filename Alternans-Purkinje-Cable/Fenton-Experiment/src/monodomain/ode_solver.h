#ifndef _ODE_SOLVER_H_
#define _ODE_SOLVER_H_

#include <cstdbool>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include "config/stim_config_hash.h"
#include "config/config_parser.h"

#include "../models_library/model_common.h"

struct ode_solver 
{

    void *handle;

    double dt;

    bool gpu;
    int gpu_id;

    double *sv;
    struct cell_model_data model_data;

    size_t pitch;

    //User provided functions
    get_cell_model_data_fn *get_cell_model_data;
    set_ode_initial_conditions_cpu_fn *set_ode_initial_conditions_cpu;
    set_ode_initial_conditions_gpu_fn *set_ode_initial_conditions_gpu;
    solve_model_ode_cpu_fn *solve_model_ode_cpu;
    solve_model_ode_gpu_fn *solve_model_ode_gpu;


};

struct ode_solver* new_ode_solver();

#endif