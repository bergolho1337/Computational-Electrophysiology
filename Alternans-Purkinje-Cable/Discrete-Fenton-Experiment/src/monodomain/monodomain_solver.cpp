#include "monodomain_solver.h"

struct monodomain_solver* new_monodomain_solver ()
{
    struct monodomain_solver *result = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));
    
    // Default values
    result->cm = 1.2;
    result->beta = 0.14;
    
    return result;
}

void configure_monodomain_from_options (struct monodomain_solver *solver, struct user_options *options)
{
    assert(solver);
    assert(options);

    double dt = options->dt;
    solver->dt = dt;

    double tmax = options->tmax;
    solver->tmax = tmax;

    int use_steady_state = options->use_steady_state;
    if (use_steady_state)
        solver->use_steady_state = true;
    else
        solver->use_steady_state = false;

    // Calculate the beta parameter based on diameter of the Purkinje cell
    double diameter = options->diameter;
    solver->beta = 4.0f / diameter * 1.0e-04;

}