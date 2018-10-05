#ifndef _ODE_SOLVER_H_
#define _ODE_SOLVER_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>

struct ode_solver
{
    int num_ode_equations;
    int n_active_cells;

    double *sv;
};

struct ode_solver* new_ode_solver ();

#endif