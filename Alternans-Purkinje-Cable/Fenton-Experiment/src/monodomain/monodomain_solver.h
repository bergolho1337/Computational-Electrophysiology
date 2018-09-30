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

#endif