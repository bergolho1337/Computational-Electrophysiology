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

    float *sv;
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
void configure_ode_solver_from_options(struct ode_solver *solver, struct user_options *options);
void init_ode_solver_with_cell_model(struct ode_solver* solver);
void set_ode_initial_conditions_for_all_volumes(struct ode_solver *solver, uint32_t num_cells);
void set_ode_initial_conditions_using_steady_state(struct ode_solver *the_ode_solver, const uint32_t n_active, char *input_steady_filename);
void solve_all_volumes_odes(struct ode_solver *the_ode_solver, uint32_t n_active, double cur_time, int num_steps,
                            struct stim_config_hash *stim_configs);

#endif