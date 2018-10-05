#ifndef _MONODOMAIN_SOLVER_H_
#define _MONODOMAIN_SOLVER_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

struct monodomain_solver
{
    double cm;
    double beta;

    double dt;
    double tmax;
    int use_steady_state;

};

struct monodomain_solver* new_monodomain_solver ();

#endif