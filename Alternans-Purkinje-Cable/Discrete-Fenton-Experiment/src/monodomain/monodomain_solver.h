#ifndef _MONODOMAIN_SOLVER_H_
#define _MONODOMAIN_SOLVER_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "../config/user_config.h"

using namespace std;

struct monodomain_solver
{
    double cm;
    double beta;

    double dt;
    double tmax;
    bool use_steady_state;

};

struct monodomain_solver* new_monodomain_solver ();
void configure_monodomain_from_options (struct monodomain_solver *solver, struct user_options *options);

#endif