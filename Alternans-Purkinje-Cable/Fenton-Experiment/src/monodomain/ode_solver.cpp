#include "ode_solver.h"
#include <cstring>
#include <dlfcn.h>
#include <cassert>
#include "../utils/logfile_utils.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

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

void set_ode_initial_conditions_for_all_volumes(struct ode_solver *solver, uint32_t num_cells) 
{

    bool get_initial_v = !std::isfinite(solver->model_data.initial_v);
    bool get_neq = solver->model_data.number_of_ode_equations == -1;

    (*(solver->get_cell_model_data))(&(solver->model_data), get_initial_v, get_neq);
    int n_odes = solver->model_data.number_of_ode_equations;

    if (solver->gpu) 
    {
#ifdef COMPILE_CUDA

        set_ode_initial_conditions_gpu_fn *soicg_fn_pt = solver->set_ode_initial_conditions_gpu;

        if(!soicg_fn_pt) 
        {
            fprintf(stderr, "The ode solver was set to use the GPU, \n "
                    "but no function called set_model_initial_conditions_gpu "
                    "was provided in the %s shared library file\n", solver->model_data.model_library_path);
            exit(11);
        }

        if(solver->sv != NULL) 
        {
            check_cuda_errors(cudaFree(solver->sv));
        }

        solver->pitch = soicg_fn_pt(&(solver->sv), num_cells);
#endif
    } 
    else 
    {

        set_ode_initial_conditions_cpu_fn *soicc_fn_pt = solver->set_ode_initial_conditions_cpu;

        if(!soicc_fn_pt) 
        {
            fprintf(stderr, "The ode solver was set to use the CPU, \n "
                    "but no function called set_model_initial_conditions_cpu "
                    "was provided in the %s shared library file\n", solver->model_data.model_library_path);
            exit(11);
        }

        if(solver->sv != NULL) 
        {
            free(solver->sv);
        }

        solver->sv = (float*)malloc(n_odes*num_cells*sizeof(float));

		int i;

        #pragma omp parallel for
        for(i = 0; i < num_cells; i++) 
        {
            soicc_fn_pt(solver->sv + (i*n_odes));
        }
    }
}

void set_ode_initial_conditions_using_steady_state(struct ode_solver *the_ode_solver, const uint32_t n_active, char *input_steady_filename)
{
    print_to_stdout_and_file("Using a Steady State solution for the initial conditions ...\n");

    int n_odes = the_ode_solver->model_data.number_of_ode_equations;
    float *sv = the_ode_solver->sv;
    size_t pitch = the_ode_solver->pitch;
    bool use_gpu = the_ode_solver->gpu;
    size_t mem_size = n_active * n_odes * sizeof (float);

    // Reading steady state solution
    FILE *file = fopen(input_steady_filename,"r");
    float *sv_sst = (float*)malloc(mem_size);
    for (int i = 0; i < n_active; i++)
        for (int j = 0; j < n_odes; j++)
            fscanf(file,"%f",&sv_sst[i*n_odes+j]);
    fclose(file);

    // Moving the steady-state solution to the CPU || GPU state-vector ...
    if (use_gpu) 
    {
#ifdef COMPILE_CUDA

        check_cuda_errors(cudaMemcpy2D(sv, n_active * sizeof(float), sv_sst, pitch, n_active * sizeof(float), n_odes, cudaMemcpyHostToDevice));
        
#endif
    } 
    else 
    {
		for (int i = 0; i < n_active; i++) 
            for (int j = 0; j < n_odes; j++)
                sv[i*n_odes+j] = sv_sst[i*n_odes+j];
    }

    free(sv_sst);

}

void solve_all_volumes_odes(struct ode_solver *the_ode_solver, uint32_t n_active, double cur_time, int num_steps,
                            struct stim_config_hash *stim_configs) 
{

    assert(the_ode_solver->sv);

    float dt = the_ode_solver->dt;
    //int n_odes = the_ode_solver->model_data.number_of_ode_equations;
    float *sv = the_ode_solver->sv;

    double time = cur_time;

    float *merged_stims = (float*)calloc(sizeof(float), n_active);

    struct stim_config *tmp = NULL;
    float stim_period;
    float stim_start, stim_duration;
    float start_period, end_period, period_step;
    int n_cycles;

    float new_time;

	int i;

    if(stim_configs) 
    {
        for (int k = 0; k < stim_configs->size; k++) 
        {
            for (struct stim_config_elt *e = stim_configs->table[k % stim_configs->size]; e != 0; e = e->next) 
            {
                tmp = e->value;
                stim_start = tmp->stim_start;
                stim_duration = tmp->stim_duration;
                start_period = tmp->start_period;
                end_period = tmp->end_period;
                period_step = tmp->period_step;
                n_cycles = tmp->n_cycles;
                
                for (int j = 0; j < num_steps; ++j) 
                {
                    new_time = 0.0f;
                    // New Jhonny stimulus protocol for alternans simulations ...
                    for (double new_period = start_period; new_period >= end_period; new_period -= period_step)
                    {
                        if ( time >= new_time && (time < new_time + n_cycles*new_period || new_period == end_period) )
                        {
                            stim_period = new_period;
                            time -= new_time;
                            break;
                        }
                        new_time += n_cycles*new_period;

                    }
                    if( (time-floor(time/stim_period)*stim_period>=stim_start) && ( time - floor(time/stim_period)*stim_period <= stim_start + stim_duration ) )
                    {
                        #pragma omp parallel for
                        for (i = 0; i < n_active; i++) 
                        {
                            merged_stims[i] = tmp->spatial_stim_currents[i];
                        }
                    }

                    // Old Sachetto's stimulus protocol ...
                    /*
                    if ((time >= stim_start) && (time <= stim_start + stim_duration)) 
                    {
                        #pragma omp parallel for
                        for (i = 0; i < n_active; i++) 
                        {
                            merged_stims[i] = tmp->spatial_stim_currents[i];
                        }
                    }
                    */
                    time += dt;
                }
                time = cur_time;
            }
        }
    }


    if(the_ode_solver->gpu) 
    {
#ifdef COMPILE_CUDA
        solve_model_ode_gpu_fn *solve_odes_pt = the_ode_solver->solve_model_ode_gpu;
        solve_odes_pt(dt, sv, merged_stims, NULL, n_active, num_steps, NULL,
                      0);

#endif
    }
    else 
    {
        solve_model_ode_cpu_fn *solve_odes_pt = the_ode_solver->solve_model_ode_cpu;
        solve_odes_pt(dt, sv, merged_stims, NULL, n_active, num_steps, NULL);
    }

    free(merged_stims);
    
}