#ifndef _MONODOMAIN_SOLVER_H_
#define _MONODOMAIN_SOLVER_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <Eigen/Sparse>
#include <vector>
#include <algorithm>
#include "../config/user_config.h"
#include "../purkinje/purkinje.h"
#include "ode_solver.h"

using namespace std;

struct monodomain_solver
{
    double cm;
    double beta;
    double sigma_c;
    double G_gap;

    double dt;
    double tmax;
    bool use_steady_state;

};

struct monodomain_solver* new_monodomain_solver ();
void configure_monodomain_from_options (struct monodomain_solver *solver, struct user_options *options);
void solve_monodomain (struct monodomain_solver *the_monodomain_solver,\
                        struct ode_solver *the_ode_solver,\
                        struct graph *the_purkinje_network,\
                        struct user_options *options);
Eigen::SparseMatrix<double> assemble_matrix (struct monodomain_solver *the_monodomain_solver,\
                      struct graph *the_purkinje_network);

void print_solver_info (struct monodomain_solver *the_monodomain_solver,\
                        struct ode_solver *the_ode_solver,\
                        struct graph *the_purkinje_network,\
                        struct user_options *options);
void write_to_VTK (struct graph *the_purkinje_network, struct ode_solver *the_ode_solver, int iter);
void assemble_matrix (struct monodomain_solver *the_monodomain_solver,\
                      struct graph *the_purkinje_network, Eigen::SparseMatrix<double> &A);
void assemble_load_vector (struct ode_solver *the_ode_solver, const double A, Eigen::VectorXd &b);
void update_monodomain (struct ode_solver *the_ode_solver, const Eigen::VectorXd vm);
void solve_all_volumes_odes (struct ode_solver *the_ode_solver, double *vstar, const double dt, const double cur_time);
void set_stimulus (double *merged_stim, const uint32_t n_cells, const uint32_t cur_time);
void next_timestep (struct ode_solver *the_ode_solver);                        

void swap (double **a, double **b);

#endif