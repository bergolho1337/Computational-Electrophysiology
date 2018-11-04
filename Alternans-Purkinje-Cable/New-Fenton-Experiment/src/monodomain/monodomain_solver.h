//
// Created by sachetto on 29/09/17.
//

#ifndef FENTON_SOLVER_H
#define FENTON_SOLVER_H

//#include "../alg/grid/grid.h"
//#include "ode_solver.h"
//#include "config/stim_config_hash.h"
//#include "config/extra_data_config.h"
#include "config/config_parser.h"

#include <cstdbool>
#include <cstdint>
#include <dlfcn.h>

#include "../models_library/model_common.h"

// Structure for the ploting of the volumes
struct plot
{
  FILE **plotFile;                            // Reference to the file of the plot volume
  int np;                                     // Number of plot volumes
  int *ids;                                   // Identifier of the plot volume
};

struct monodomain_solver 
{

    bool use_steady_state;
    int num_threads;
    double final_time;

    double beta, cm; // micrometers
    double sigma_c;
    double G_gap;

    double dt;

    double start_h;
    double start_diameter;


    struct plot *plot;

    void *handle;
    struct cell_model_data model_data;
    
    //User provided functions
    get_cell_model_data_fn *get_cell_model_data;
    set_ode_initial_conditions_cpu_fn *set_ode_initial_conditions_cpu;
    //set_ode_initial_conditions_gpu_fn *set_ode_initial_conditions_gpu;
    solve_model_ode_cpu_fn *solve_model_ode_cpu;
    //solve_model_ode_gpu_fn *solve_model_ode_gpu;
};

struct monodomain_solver *new_monodomain_solver ();

void configure_monodomain_solver_from_options(struct monodomain_solver *the_monodomain_solver,
                                              struct user_options *options);
void configure_plot_cells (struct monodomain_solver *the_monodomain_solver,
                                              struct user_options *options);

void solve_monodomain(struct monodomain_solver *the_monodomain_solver,
                      struct grid *the_grid, struct user_options *configs);

void set_celular_model (struct monodomain_solver *solver, 
                        struct user_options *configs);
/*
void save_old_cell_positions (struct grid *the_grid);
void update_cells_to_solve (struct grid *the_grid, struct ode_solver *solver);
void set_initial_conditions_all_volumes (struct monodomain_solver *the_solver, struct grid *the_grid, double initial_v);

void print_solver_info(struct monodomain_solver *the_monodomain_solver, struct ode_solver *the_ode_solver,
                       struct grid *the_grid, struct user_options *options);

void update_ode_state_vector(struct ode_solver *the_ode_solver, struct grid *the_grid, uint32_t max_number_of_cells);

void set_ode_extra_data(struct extra_data_config *config, struct grid *the_grid, struct ode_solver *the_ode_solver);
void set_spatial_stim(struct stim_config_hash *stim_configs, struct grid *the_grid);

void update_monodomain(uint32_t initial_number_of_cells, uint32_t num_active_cells, struct cell_node **active_cells,
                       double beta,
                       double cm, double dt_edp, real *sv, int n_equations_cell_model, bool use_gpu);




bool print_result(const struct grid *the_grid, const struct user_options *configs, int count, bool save_in_binary);

*/
#endif // MONOALG3D_SOLVER_H
