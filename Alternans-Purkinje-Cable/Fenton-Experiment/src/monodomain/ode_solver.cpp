#include "ode_solver.h"

struct ode_solver* new_ode_solver()
{
    struct ode_solver* result = (struct ode_solver *) malloc(sizeof(struct ode_solver));
    result->sv = NULL;
    result->handle = NULL;

    //result->get_cell_model_data = NULL;
    //result->set_ode_initial_conditions_cpu = NULL;
    //result->solve_model_ode_cpu = NULL;

    //result->set_ode_initial_conditions_gpu = NULL;
    //result->solve_model_ode_gpu = NULL;
    //result->model_data.initial_v = INFINITY;
    //result->model_data.number_of_ode_equations = -1;

    return result;
}