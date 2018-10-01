#ifndef _MONODOMAIN_SOLVER_H_
#define _MONODOMAIN_SOLVER_H_

#include <cstdbool>
#include <cstdint>

#include "config/config_parser.h"

#include <Eigen/Sparse>

struct monodomain_solver 
{

    int num_threads;
    double final_time;

    double beta, cm; // micrometers
    double sigma;

    // Time used for solving wave equation.
    double dt;

    // TO DO: Put Eigen matrix here ...


};

struct monodomain_solver *new_monodomain_solver ();
void configure_monodomain_solver_from_options (struct monodomain_solver *the_monodomain_solver,
                                               struct user_options *options);
                                               
void solve_monodomain (struct monodomain_solver *monodomain_solver, struct ode_solver *ode_solver,\
                        struct grid *grid, struct user_options *configs);

void set_spatial_purkinje (struct purkinje_config *pk_config, struct grid *grid);

void print_solver_info (struct monodomain_solver *the_monodomain_solver, struct ode_solver *the_ode_solver,
                        struct grid *the_grid, struct user_options *options);

void set_initial_conditions_all_volumes (struct monodomain_solver *the_solver, struct grid *the_grid, double initial_v);

void set_spatial_stim(struct stim_config_hash *stim_configs, struct grid *the_grid);
void set_stimulus (struct stim_config *tmp, struct grid *grid);

void update_ode_state_vector (struct ode_solver *the_ode_solver, struct grid *the_grid, uint32_t max_number_of_cells);

void update_monodomain (uint32_t initial_number_of_cells, uint32_t num_active_cells, struct cell_node **active_cells,
                        double beta, double cm, double dt_edp, float *sv, int n_equations_cell_model, bool use_gpu);

bool print_result(const struct grid *the_grid, const struct user_options *configs, int count, bool save_in_binary);

Eigen::SparseMatrix<double> assembly_matrix (struct monodomain_solver *monodomain_solver, struct grid *grid, struct purkinje_config *pk_config);
void assembly_load_vector (Eigen::VectorXd &b, struct cell_node **active_cells, uint32_t num_active_cells);

void move_solution_to_cells (Eigen::VectorXd &x, struct cell_node **active_cells, uint32_t num_active_cells);

#endif