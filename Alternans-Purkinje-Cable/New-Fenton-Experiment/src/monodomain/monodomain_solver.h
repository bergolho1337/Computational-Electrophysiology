//
// Created by sachetto on 29/09/17.
//

#ifndef FENTON_SOLVER_H
#define FENTON_SOLVER_H

//#include "../alg/grid/grid.h"
//#include "ode_solver.h"
//#include "config/extra_data_config.h"
#include "config/stim_config_hash.h"
#include "config/config_parser.h"
#include "config/config_helper.h"

#include <cstdbool>
#include <cstdint>
#include <dlfcn.h>

#include <Eigen/Sparse>
#include <vector>

#include "../models_library/model_common.h"

// Structure for the ploting of the volumes
struct plot
{
  FILE **plotFile;                            // Reference to the file of the plot volume
  int np;                                     // Number of plot volumes
  int *ids;                                   // Identifier of the plot volume
};

struct derivative
{
  double t;                                   // Time of the maximum derivative
  double value;                               // Derivative value
};

struct velocity
{
  FILE *velocity_file;                         // Reference to the file where the velocity will be stored
  int np;                                     // Number of volumes that the velocity will be calculated
  int id_source;                              // Identifier of the source volume
  int *ids;                                   // Identifier of the sink points (ids[0] -> source) 
  double t1;                                  // Time when the despolarization occur on the source volume 
  double *t2;                                 // Time when the despolarization occur on the sink volumes
};


struct monodomain_solver 
{

    bool use_steady_state;
    int num_threads;

    double beta, cm; // micrometers
    double sigma_c;
    double G_gap;

    double dt;
    double final_time;
    int M;

    double start_h;
    double start_diameter;

    int num_volumes;

    struct plot *plot;
    struct derivative *dvdt;
    struct velocity *vel;
    struct control_volume *volumes;

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

void set_initial_conditions_from_file (struct monodomain_solver *solver, struct user_options *options);
void set_initial_conditions_default (struct monodomain_solver *solver);

void set_matrix (Eigen::SparseMatrix<double> &a,\
                struct monodomain_solver *solver,\
                struct grid *the_grid);
void assemble_load_vector (Eigen::VectorXd &b,\
                            struct monodomain_solver *solver,\
                            struct grid *the_grid);
void move_v_star (const Eigen::VectorXd vm,\
                    struct monodomain_solver *solver);

void set_spatial_stim(struct stim_config_hash *stim_configs,\
                      struct grid *the_grid);

void solve_odes (const double t,\
                 struct monodomain_solver *solver,
                 struct stim_config_hash *stim_configs);

double* merge_stimulus (struct stim_config_hash *stim_configs,\
                    const int np, const double cur_time);

void calc_max_derivative (struct monodomain_solver *solver,\
                        const double t, const double start_period);

void calc_velocity (struct monodomain_solver *solver, struct grid *the_grid);

void next_timestep (struct monodomain_solver *solver);
void swap (double **a, double **b);

void print_solver_info (struct monodomain_solver *the_monodomain_solver,\
                        struct grid *the_grid,\
                        struct user_options *configs);
void print_progress (int iter, int max_iter);

void write_plot_data (struct monodomain_solver *solver, double t);
void write_steady_state_to_file (FILE *sst_file,\
                                struct monodomain_solver *solver);
void write_VTK_to_file (struct monodomain_solver *solver,\
                        struct grid *the_grid, int iter);

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
