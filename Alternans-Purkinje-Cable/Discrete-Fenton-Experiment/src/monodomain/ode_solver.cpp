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

    // Using Noble 1962 celular model
    strcpy(the_ode_solver->model_name,"Noble 1962");

    // Get the number of equations from the model
    get_cell_model_data_fn *cell_model_data = the_ode_solver->get_cell_model_data;
    cell_model_data = init_cell_model_data;

    cell_model_data(&the_ode_solver->num_ode_equations);
    the_ode_solver->n_active_cells = num_volumes;

    // Get the initial condition function
    the_ode_solver->set_ode_initial_conditions_cpu = set_model_initial_conditions_cpu;

    // Get the ODE solution function
    the_ode_solver->solve_model_ode_cpu = solve_model_odes_cpu;

    // Allocate memory
    uint32_t n_odes = the_ode_solver->num_ode_equations;
    the_ode_solver->sv = (double*)malloc(sizeof(double)*num_volumes*n_odes);

}

void set_ode_initial_condition_for_all_volumes (struct ode_solver *the_ode_solver)
{
    // TO DO: Include the steady-state initial condition

    uint32_t n_odes = the_ode_solver->num_ode_equations;
    uint32_t n_active = the_ode_solver->n_active_cells;
    set_ode_initial_conditions_cpu_fn *soicf_fn_pt = the_ode_solver->set_ode_initial_conditions_cpu;

    uint32_t i;

    #pragma omp parallel for
    for(i = 0; i < n_active; i++) 
    {
        soicf_fn_pt(the_ode_solver->sv + (i*n_odes));
    }
}

void print_state_vector (struct ode_solver *the_ode_solver)
{
    uint32_t n_odes = the_ode_solver->num_ode_equations;
    uint32_t n_active = the_ode_solver->n_active_cells;
    double *sv = the_ode_solver->sv;

    for (uint32_t i = 0; i < n_active; i++)
    {
        for (uint32_t j = 0; j < n_odes-1; j++)
        {
            printf("%.10lf ",sv[i*n_odes+j]);
        }
        printf("%.10lf\n",sv[i*n_odes+(n_odes-1)]);
    }
}