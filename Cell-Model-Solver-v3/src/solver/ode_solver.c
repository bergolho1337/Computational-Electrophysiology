//
// Created by sachetto on 02/10/17.
//

#include "ode_solver.h"

#include <string.h>
#ifdef _MSC_VER
#include "../dlfcn-win32/dlfcn.h"
#else
#include <dlfcn.h>
#endif
#include <assert.h>
#include "../utils/logfile_utils.h"



struct ode_solver* new_ode_solver() 
{
    struct ode_solver* result = (struct ode_solver *) malloc(sizeof(struct ode_solver));
    result->sv = NULL;
    result->cells_to_solve = NULL;
    result->handle = NULL;

    result->get_cell_model_data = NULL;
    result->set_ode_initial_conditions_cpu = NULL;
    result->solve_model_ode_cpu = NULL;

    result->model_data.initial_v = INFINITY;
    result->model_data.number_of_ode_equations = -1;

    result->edo_extra_data = NULL;
    result->edo_extra_data = 0;

    return result;
}

void free_ode_solver(struct ode_solver *solver) 
{
    if(solver->sv) 
    {
        free(solver->sv);
    }

    if(solver->edo_extra_data) 
    {
        free(solver->edo_extra_data);
    }

    if(solver->cells_to_solve) 
    {
        free(solver->cells_to_solve);
    }

    if(solver->model_data.model_library_path) 
    {
        free(solver->model_data.model_library_path);
    }

    if(solver->handle) 
    {
        dlclose(solver->handle);
    }

    free(solver);

}

void configure_ode_solver_from_options(struct ode_solver *solver, struct user_options *options) {
    solver->gpu_id = options->gpu_id;
    solver->min_dt = (real)options->dt_edo;
    solver->gpu = options->gpu;

    if(options->model_file_path) {
        solver->model_data.model_library_path = strdup(options->model_file_path);
    }

}

void solve_celular_model (struct ode_solver *solver, struct user_options *options)
{
    assert (options);
    assert (solver);

    print_to_stdout_and_file (LOG_LINE_SEPARATOR);

    long ode_total_time = 0, total_write_time = 0, total_config_time = 0;

    struct stop_watch solver_time, ode_time, write_time, config_time;

    init_stop_watch (&config_time);

    start_stop_watch (&config_time);

    // _____________________________________________________________
    // ====== Begin configuation =============
    init_ode_solver_with_cell_model (solver);
    struct stim_config_hash *stimuli_configs = options->stim_configs;

    double last_stimulus_time = -1.0;

    if (stimuli_configs) 
    {
        // Init all stimuli
        STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE_AND_KEY (stimuli_configs, init_stim_functions);

        //Find last stimuli
        size_t s_size = stimuli_configs->size;
        double s_end;
        for (int i = 0; i < s_size; i++) 
        {
            for (struct stim_config_elt *e = stimuli_configs->table[i % s_size]; e != 0; e = e->next) 
            {
                s_end = e->value->stim_start + e->value->stim_duration;
                if(s_end > last_stimulus_time) last_stimulus_time = s_end;
            }
        }
    }
    // ====== End configuation =============
    
    int count = 0;
    bool save_to_file = (options->out_dir_name != NULL);
    double dt_edo = solver->min_dt;

    print_to_stdout_and_file ("Setting ODE's initial conditions\n");
    set_ode_initial_conditions(solver);


}

void init_ode_solver_with_cell_model(struct ode_solver* solver) {

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

    solver->get_cell_model_data = dlsym(solver->handle, "init_cell_model_data");
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "init_cell_model_data function not found in the provided model library\n");
        if(!isfinite(solver->model_data.initial_v)) {
            fprintf(stderr, "initial_v not provided in the [cell_model] of the config file! Exiting\n");
            exit(1);
        }

    }

    solver->set_ode_initial_conditions_cpu = dlsym(solver->handle, "set_model_initial_conditions_cpu");
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "set_model_initial_conditions function not found in the provided model library\n");
        exit(1);
    }

    solver->solve_model_ode_cpu = dlsym(solver->handle, "solve_model_odes_cpu");
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "solve_model_odes_cpu function not found in the provided model library\n");
        exit(1);
    }

}

void set_ode_initial_conditions (struct ode_solver *solver) 
{

    bool get_initial_v = !isfinite(solver->model_data.initial_v);
    bool get_neq = solver->model_data.number_of_ode_equations == -1;

    (*(solver->get_cell_model_data))(&(solver->model_data), get_initial_v, get_neq);
    int n_odes = solver->model_data.number_of_ode_equations;
 
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

    solver->sv = (real*)malloc(n_odes*sizeof(real));

    soicc_fn_pt(solver->sv);

}

/*
void solve_all_volumes_odes(struct ode_solver *the_ode_solver, uint32_t n_active, double cur_time, int num_steps,
                            struct stim_config_hash *stim_configs) {

    assert(the_ode_solver->sv);

    real dt = the_ode_solver->min_dt;
    //int n_odes = the_ode_solver->model_data.number_of_ode_equations;
    real *sv = the_ode_solver->sv;

    void *extra_data = the_ode_solver->edo_extra_data;
    size_t extra_data_size = the_ode_solver->extra_data_size;

    double time = cur_time;

    real *merged_stims = (real*)calloc(sizeof(real), n_active);

    struct stim_config *tmp = NULL;
    real stim_start, stim_dur;

	int i;

    if(stim_configs) {
        for (int k = 0; k < stim_configs->size; k++) {
            for (struct stim_config_elt *e = stim_configs->table[k % stim_configs->size]; e != 0; e = e->next) {
                tmp = e->value;
                stim_start = tmp->stim_start;
                stim_dur = tmp->stim_duration;
                for (int j = 0; j < num_steps; ++j) {
                    if ((time >= stim_start) && (time <= stim_start + stim_dur)) {
                        #pragma omp parallel for
                        for (i = 0; i < n_active; i++) {
                            merged_stims[i] = tmp->spatial_stim_currents[i];
                        }
                    }
                    time += dt;
                }
                time = cur_time;
            }
        }
    }


    if(the_ode_solver->gpu) {
#ifdef COMPILE_CUDA
        solve_model_ode_gpu_fn *solve_odes_pt = the_ode_solver->solve_model_ode_gpu;
        solve_odes_pt(dt, sv, merged_stims, the_ode_solver->cells_to_solve, n_active, num_steps, extra_data,
                      extra_data_size);

#endif
    }
    else {
        solve_model_ode_cpu_fn *solve_odes_pt = the_ode_solver->solve_model_ode_cpu;
        solve_odes_pt(dt, sv, merged_stims, the_ode_solver->cells_to_solve, n_active, num_steps, extra_data);
    }

    free(merged_stims);
}

void update_state_vectors_after_refinement(struct ode_solver *ode_solver, const uint32_t *refined_this_step) {

    assert(ode_solver);
    assert(ode_solver->sv);

    size_t num_refined_cells = sb_count(refined_this_step)/8;

    real *sv = ode_solver->sv;
    int neq = ode_solver->model_data.number_of_ode_equations;
    real *sv_src;
    real *sv_dst;
    size_t  i;


    if(ode_solver->gpu) {
#ifdef COMPILE_CUDA

        size_t pitch_h = ode_solver->pitch;

		#pragma omp parallel for private(sv_src, sv_dst)
        for (i = 0; i < num_refined_cells; i++) {

            size_t index_id = i * 8;

            uint32_t index = refined_this_step[index_id];
            sv_src = &sv[index];

            for (int j = 1; j < 8; j++) {
                index = refined_this_step[index_id + j];
                sv_dst = &sv[index];
                check_cuda_errors(cudaMemcpy2D(sv_dst, pitch_h, sv_src, pitch_h, sizeof(real), (size_t )neq, cudaMemcpyDeviceToDevice));
            }


        }
        //TODO: test if is faster to update the GPU using a kernel or a host function with cudaMemcpy2D
        //ode_solver->update_gpu_fn(sv, refined_this_step->base, num_refined_cells, neq);

#endif
    }
    else {

        #pragma omp parallel for private(sv_src, sv_dst)
        for (i = 0; i < num_refined_cells; i++) {

            size_t index_id = i * 8;

            uint32_t index = refined_this_step[index_id];
            sv_src = &sv[index * neq];

            for (int j = 1; j < 8; j++) {
                index = refined_this_step[index_id + j];
                sv_dst = &sv[index * neq];
                memcpy(sv_dst, sv_src, neq * sizeof(real));
            }


        }
    }

}
*/