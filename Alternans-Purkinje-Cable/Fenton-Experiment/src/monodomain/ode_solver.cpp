#include "ode_solver.h"
#include <cstring>
#include <dlfcn.h>
#include <cassert>
#include "../utils/logfile_utils.h"

struct ode_solver* new_ode_solver()
{
    struct ode_solver* result = (struct ode_solver *) malloc(sizeof(struct ode_solver));
    result->sv = NULL;
    result->handle = NULL;

    result->get_cell_model_data = NULL;
    result->set_ode_initial_conditions_cpu = NULL;
    result->solve_model_ode_cpu = NULL;

    result->set_ode_initial_conditions_gpu = NULL;
    result->solve_model_ode_gpu = NULL;
    result->model_data.initial_v = INFINITY;
    result->model_data.number_of_ode_equations = -1;

    return result;
}

void configure_ode_solver_from_options(struct ode_solver *solver, struct user_options *options) 
{
    solver->gpu_id = options->gpu_id;
    solver->dt = (float)options->dt;
    solver->gpu = options->gpu;

    if(options->model_file_path) 
    {
        solver->model_data.model_library_path = strdup(options->model_file_path);
    }

}

void init_ode_solver_with_cell_model(struct ode_solver* solver) 
{

    char *error;

    if(!solver->model_data.model_library_path) {
        fprintf(stderr, "model_library_path not provided. Exiting!\n");
        exit(1);
    }

    print_to_stdout_and_file("Opening %s as model lib\n", solver->model_data.model_library_path);

    solver->handle = dlopen (solver->model_data.model_library_path, RTLD_LAZY);
    if (!solver->handle) {
        fprintf(stderr, "%s\n", dlerror());
        exit(1);
    }

    solver->get_cell_model_data = (get_cell_model_data_fn*)dlsym(solver->handle, "init_cell_model_data");
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "init_cell_model_data function not found in the provided model library\n");
        if(!std::isfinite(solver->model_data.initial_v)) 
        {
            fprintf(stderr, "intial_v not provided in the [cell_model] of the config file! Exiting\n");
            exit(1);
        }

    }

    solver->set_ode_initial_conditions_cpu = (set_ode_initial_conditions_cpu_fn*)dlsym(solver->handle, "set_model_initial_conditions_cpu");
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "set_model_initial_conditions function not found in the provided model library\n");
        exit(1);
    }

    solver->solve_model_ode_cpu = (solve_model_ode_cpu_fn*)dlsym(solver->handle, "solve_model_odes_cpu");
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "solve_model_odes_cpu function not found in the provided model library\n");
        exit(1);
    }

#ifdef COMPILE_CUDA
    solver->set_ode_initial_conditions_gpu = (set_ode_initial_conditions_gpu_fn*)dlsym(solver->handle, "set_model_initial_conditions_gpu");
    if ((error = dlerror()) != NULL)  {
        fputs(error, stderr);
        fprintf(stderr, "set_model_initial_conditions_gpu function not found in the provided model library\n");
        exit(1);
    }

    solver->solve_model_ode_gpu = (solve_model_ode_gpu_fn*)dlsym(solver->handle, "solve_model_odes_gpu");
    if ((error = dlerror()) != NULL)  {
        fputs(error, stderr);
        fprintf(stderr, "\nsolve_model_odes_gpu function not found in the provided model library\n");
        exit(1);
    }
    
#endif

}