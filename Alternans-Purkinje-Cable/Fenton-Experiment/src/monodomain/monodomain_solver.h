#ifndef _MONODOMAIN_SOLVER_H_
#define _MONODOMAIN_SOLVER_H_

#include <cstdbool>
#include <cstdint>

#include "config/config_parser.h"

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

#endif