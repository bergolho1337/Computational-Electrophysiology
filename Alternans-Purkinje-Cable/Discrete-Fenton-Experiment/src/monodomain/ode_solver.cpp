#include "ode_solver.h"

struct ode_solver* new_ode_solver ()
{
    struct ode_solver* result = (struct ode_solver*)malloc(sizeof(struct ode_solver));
    
    result->n_active_cells = 0;
    result->num_ode_equations = -1;

    result->sv = NULL;
    result->get_cell_model_data = NULL;
    result->set_ode_initial_conditions_cpu = NULL;
    result->solve_model_ode_cpu = NULL;

    return result;
}

void configure_ode_solver (struct ode_solver *the_ode_solver, const uint32_t num_volumes)
{
    assert(the_ode_solver);

    // Get the number of equations from the model
    get_cell_model_data_fn *cell_model_data = the_ode_solver->get_cell_model_data;
    cell_model_data = init_cell_model_data;

    cell_model_data(&the_ode_solver->num_ode_equations);
    the_ode_solver->n_active_cells = num_volumes;

    // Get the initial condition function
    set_ode_initial_conditions_cpu_fn *model_initial_conditions = the_ode_solver->set_ode_initial_conditions_cpu;
    model_initial_conditions = set_model_initial_conditions_cpu;

    // Get the ODE solution function
    solve_model_ode_cpu_fn *solve_model_ode = the_ode_solver->solve_model_ode_cpu;
    solve_model_ode = solve_model_odes_cpu;


}