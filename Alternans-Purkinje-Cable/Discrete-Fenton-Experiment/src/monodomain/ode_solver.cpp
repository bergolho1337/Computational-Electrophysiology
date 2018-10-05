#include "ode_solver.h"

struct ode_solver* new_ode_solver ()
{
    struct ode_solver* result = (struct ode_solver*)malloc(sizeof(struct ode_solver));
    
    result->n_active_cells = 0;
    result->num_ode_equations = -1;
    result->sv = NULL;

    return result;
}