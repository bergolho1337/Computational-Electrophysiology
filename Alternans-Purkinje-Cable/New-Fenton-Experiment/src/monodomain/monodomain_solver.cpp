//
// Created by sachetto on 03/10/17.
//

#include "monodomain_solver.h"
//#include "../utils/logfile_utils.h"
//#include "../utils/stop_watch.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cinttypes>
//#include "../string/sds.h"


//#include "config/domain_config.h"
//#include "config/purkinje_config.h"
//#include "config/assembly_matrix_config.h"
//#include "config/linear_system_solver_config.h"

/*
static inline double ALPHA (double beta, double cm, double dt, double h) 
{
    return (((beta * cm) / dt) * UM2_TO_CM2) * pow (h, 3.0);
}
*/

struct monodomain_solver *new_monodomain_solver () 
{

    struct monodomain_solver *result = (struct monodomain_solver *)malloc(sizeof (struct monodomain_solver));

    result->cm = 1.2;

    return result;
}

void configure_monodomain_solver_from_options (struct monodomain_solver *the_monodomain_solver,
                                               struct user_options *options) 
{

    assert (the_monodomain_solver);
    assert (options);

    the_monodomain_solver->num_threads = options->num_threads;
    the_monodomain_solver->final_time = options->final_time;
    the_monodomain_solver->use_steady_state = options->use_steady_state;

    the_monodomain_solver->dt = options->dt;
    the_monodomain_solver->start_h = options->start_h;
    the_monodomain_solver->start_diameter = options->start_diameter;
    the_monodomain_solver->sigma_c = options->sigma_c;
    the_monodomain_solver->G_gap = options->G_gap;
    
    the_monodomain_solver->beta = 4.0 / options->start_diameter * 1.0e-04;

    configure_plot_cells(the_monodomain_solver,options);
}

void configure_plot_cells (struct monodomain_solver *the_monodomain_solver,
                                              struct user_options *options)
{
    std::string plot_filename = options->plot_filename;
    
    // Allocate and initialize the vector 'ids' with the volumes to be plotted from the input file .plt
    the_monodomain_solver->plot = (struct plot*)malloc(sizeof(struct plot));
    FILE *pltFile = fopen(plot_filename.c_str(),"r");
    if (!pltFile)
    {
        std::cerr << "[-] ERROR! Cannot open PLT file!" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!fscanf(pltFile,"%d",&the_monodomain_solver->plot->np)) 
    {
        std::cerr << "[-] ERROR! Reading PLT file!" << std::endl;
        exit(EXIT_FAILURE);
    }
    the_monodomain_solver->plot->ids = (int*)malloc(sizeof(int)*the_monodomain_solver->plot->np);

    for (int i = 0; i < the_monodomain_solver->plot->np; i++)
    {
        if (!fscanf(pltFile,"%d",&the_monodomain_solver->plot->ids[i]))
        {
            std::cerr << "[-] ERROR! Reading PLT file!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    fclose(pltFile);
}

void solve_monodomain(struct monodomain_solver *the_monodomain_solver,
                      struct grid *the_grid, struct user_options *configs)
{
    assert(the_monodomain_solver);
    assert(the_grid);
    assert(configs);

    // Main configuration
    set_celular_model(the_monodomain_solver,configs);
}

void set_celular_model (struct monodomain_solver *solver, 
                        struct user_options *configs)
{
    void *handle = dlopen (configs->model_file_path.c_str(), RTLD_LAZY);
    if (!handle) 
    {
        fprintf(stderr, "%s\n", dlerror());
        exit(EXIT_FAILURE);
    }
    char *error;
    solver->get_cell_model_data = (get_cell_model_data_fn*)dlsym(handle,"get_cell_model_data");
    if ((error = dlerror()) != NULL)  
    {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "'get_cell_model_data' not found in the provided model library\n");
        exit(EXIT_FAILURE);
    }
    
    solver->set_ode_initial_conditions_cpu = (set_ode_initial_conditions_cpu_fn*)dlsym(handle,"set_model_initial_conditions_cpu");
    if ((error = dlerror()) != NULL)  
    {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "'set_ode_initial_conditions_cpu' not found in the provided model library\n");
        exit(EXIT_FAILURE);
    }

    solver->solve_model_ode_cpu = (solve_model_ode_cpu_fn*)dlsym(handle,"solve_model_odes_cpu");
    if ((error = dlerror()) != NULL)  
    {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "'solve_model_ode_cpu' not found in the provided model library\n");
        exit(EXIT_FAILURE);
    }
}

/*
void solve_monodomain (struct monodomain_solver *the_monodomain_solver, struct ode_solver *the_ode_solver,
                       struct grid *the_grid, struct user_options *configs) {

    assert (configs);

    assert (the_grid);
    assert (the_monodomain_solver);
    assert (the_ode_solver);

#ifdef COMPILE_OPENGL
    if(configs->draw) {
        grid_to_draw = the_grid;
    }
#endif

    print_to_stdout_and_file (LOG_LINE_SEPARATOR);

    long ode_total_time = 0, cg_total_time = 0, total_write_time = 0, total_mat_time = 0, total_ref_time = 0,
         total_deref_time = 0, cg_partial, total_config_time = 0;

    uint32_t total_cg_it = 0;

    struct stop_watch solver_time, ode_time, cg_time, part_solver, part_mat, write_time, ref_time, deref_time,
        config_time;

    init_stop_watch (&config_time);

    start_stop_watch (&config_time);

    ///////MAIN CONFIGURATION BEGIN//////////////////
    init_ode_solver_with_cell_model (the_ode_solver);
    struct stim_config_hash *stimuli_configs = configs->stim_configs;
    struct extra_data_config *extra_data_config = configs->extra_data_config;
    //struct domain_config *domain_config = configs->domain_config;
    struct purkinje_config *purkinje_config = configs->purkinje_config;
    struct assembly_matrix_config *assembly_matrix_config = configs->assembly_matrix_config;
    struct linear_system_solver_config *linear_system_solver_config = configs->linear_system_solver_config;
    bool has_extra_data = (extra_data_config != NULL);

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

    // Configure the functions and set the Purkinje mesh domain
    if (purkinje_config) 
    {
        init_purkinje_functions(purkinje_config);
        purkinje_config->set_spatial_purkinje(purkinje_config,the_grid);
    } 
    else 
    {
        print_to_stdout_and_file ("No Purkinje configuration provided! Exiting!\n");
        exit (EXIT_FAILURE);
    }

    if (assembly_matrix_config) 
    {
        init_assembly_matrix_functions (assembly_matrix_config);
    } 
    else 
    {
        print_to_stdout_and_file ("No assembly matrix configuration provided! Exiting!\n");
        exit (EXIT_FAILURE);
    }

    if (linear_system_solver_config) 
    {
        init_linear_system_solver_functions (linear_system_solver_config);
    } 
    else 
    {
        print_to_stdout_and_file ("No linear solver configuration provided! Exiting!\n");
        exit (EXIT_FAILURE);
    }

    if (has_extra_data) 
    {
        init_extra_data_functions (extra_data_config);
    }
    ///////MAIN CONFIGURATION END//////////////////

    //int refine_each = the_monodomain_solver->refine_each;
    //int derefine_each = the_monodomain_solver->derefine_each; 

    //bool redo_matrix;

    bool activity;

    bool gpu = the_ode_solver->gpu;

    int count = 0;

    //double refinement_bound = the_monodomain_solver->refinement_bound;
    //double derefinement_bound = the_monodomain_solver->derefinement_bound;

    //double start_h = domain_config->start_h;
    //double max_h = domain_config->max_h;

    //double start_h = purkinje_config->start_h;
    //double max_h = purkinje_config->start_h;

    bool adaptive = the_grid->adaptive;
    //double start_adpt_at = the_monodomain_solver->start_adapting_at;
    bool save_to_file = (configs->out_dir_name != NULL);

    double dt_edp = the_monodomain_solver->dt;
    double finalT = the_monodomain_solver->final_time;

    double beta = the_monodomain_solver->beta;
    double cm = the_monodomain_solver->cm;

    double dt_edo = the_ode_solver->min_dt;

    
#ifdef COMPILE_CUDA
    if (gpu) {
        int device_count;
        int device = the_ode_solver->gpu_id;
        check_cuda_errors (cudaGetDeviceCount (&device_count));
        struct cudaDeviceProp prop;
        check_cuda_errors (cudaGetDeviceProperties (&prop, the_ode_solver->gpu_id));
        print_to_stdout_and_file ("%d devices available, running on Device %d: %s\n", device_count, device, prop.name);
        check_cuda_errors (cudaSetDevice (device));
    }
#endif

    order_grid_cells (the_grid);
    uint32_t original_num_cells = the_grid->num_active_cells;

    save_old_cell_positions (the_grid);

    if (adaptive) 
    {
        update_cells_to_solve (the_grid, the_ode_solver);
    }

    print_to_stdout_and_file ("Setting ODE's initial conditions\n");
    set_ode_initial_conditions_for_all_volumes (the_ode_solver, the_grid->num_active_cells);

    double initial_v = the_ode_solver->model_data.initial_v;

    total_config_time = stop_stop_watch (&config_time);

    print_solver_info (the_monodomain_solver, the_ode_solver, the_grid, configs);

    int ode_step = 1;

    if (dt_edp >= dt_edo) 
    {
        ode_step = (int)(dt_edp / dt_edo);
        print_to_stdout_and_file ("Solving EDO %d times before solving PDE\n", ode_step);
    } 
    else 
    {
        print_to_stdout_and_file ("WARNING: EDO time step is greater than PDE time step. Adjusting to EDO time "
                                  "step: %lf\n",
                                  dt_edo);
        dt_edp = dt_edo;
    }

    fflush (stdout);

    init_stop_watch (&solver_time);
    init_stop_watch (&ode_time);
    init_stop_watch (&cg_time);
    init_stop_watch (&part_solver);
    init_stop_watch (&part_mat);
    init_stop_watch (&write_time);
    init_stop_watch (&ref_time);
    init_stop_watch (&deref_time);

    print_to_stdout_and_file ("Assembling Monodomain Matrix Begin\n");
    start_stop_watch (&part_mat);

    set_initial_conditions_all_volumes (the_monodomain_solver, the_grid, initial_v);
    
    assembly_matrix_config->assembly_matrix(assembly_matrix_config, the_monodomain_solver, the_grid);
    
    total_mat_time = stop_stop_watch (&part_mat);
    print_to_stdout_and_file ("Assembling Monodomain Matrix End\n");
    print_to_stdout_and_file (LOG_LINE_SEPARATOR);

    start_stop_watch (&solver_time);

    int print_rate = configs->print_rate;

    bool abort_on_no_activity = the_monodomain_solver->abort_on_no_activity;
    double solver_error;
    uint32_t solver_iterations = 0;

    if (stimuli_configs)
        set_spatial_stim(stimuli_configs, the_grid);

    if (has_extra_data)
        set_ode_extra_data (extra_data_config, the_grid, the_ode_solver);

    bool save_in_binary = configs->binary;

    double cur_time = 0.0;

    print_to_stdout_and_file ("Starting simulation\n");

    while (cur_time <= finalT) 
    {

#ifdef COMPILE_OPENGL
        redraw  = count % print_rate == 0; //redraw grid
#endif

        if (save_to_file) 
        {

            if (count % print_rate == 0) 
            {
                start_stop_watch (&write_time);

                activity = print_result(the_grid, configs, count, save_in_binary);


                total_write_time += stop_stop_watch (&write_time);

                if (abort_on_no_activity) 
                {
                    if (!activity) 
                    {
                        print_to_stdout_and_file ("No activity, aborting simulation\n");
                        break;
                    }
                }
            }
        }

        if (cur_time > 0.0) 
        {
            update_ode_state_vector (the_ode_solver, the_grid, original_num_cells);
        }

        start_stop_watch (&ode_time);

        solve_all_volumes_odes (the_ode_solver, the_grid->num_active_cells, cur_time, ode_step, stimuli_configs);

        update_monodomain (original_num_cells, the_grid->num_active_cells, the_grid->active_cells, beta, cm, dt_edp,
                           the_ode_solver->sv, the_ode_solver->model_data.number_of_ode_equations, gpu);

        ode_total_time += stop_stop_watch (&ode_time);

        start_stop_watch (&cg_time);

        linear_system_solver_config->solve_linear_system(linear_system_solver_config, the_grid, &solver_iterations, &solver_error);

        cg_partial = stop_stop_watch (&cg_time);

        cg_total_time += cg_partial;

        total_cg_it += solver_iterations;

        if (count % print_rate == 0) 
        {
            print_to_stdout_and_file ("t = %lf, Iterations = "
                                      "%" PRIu32 ", Error Norm = %e, Number of Cells:"
                                      "%" PRIu32 ", Iterations time: %ld us\n",
                                      cur_time, solver_iterations, solver_error, the_grid->num_active_cells, cg_partial);
        }

        count++;
        cur_time += dt_edp;

    }

    print_to_stdout_and_file ("Resolution Time: %ld μs\n", stop_stop_watch (&solver_time));
    print_to_stdout_and_file ("ODE Total Time: %ld μs\n", ode_total_time);
    print_to_stdout_and_file ("CG Total Time: %ld μs\n", cg_total_time);
    print_to_stdout_and_file ("Mat time: %ld μs\n", total_mat_time);
    print_to_stdout_and_file ("Refine time: %ld μs\n", total_ref_time);
    print_to_stdout_and_file ("Derefine time: %ld μs\n", total_deref_time);
    print_to_stdout_and_file ("Write time: %ld μs\n", total_write_time);
    print_to_stdout_and_file ("Initial configuration time: %ld μs\n", total_config_time);
    print_to_stdout_and_file ("CG Total Iterations: %u\n", total_cg_it);
}

bool print_result(const struct grid *the_grid, const struct user_options *configs, int count, bool save_in_binary) {
    bool activity;
    sds tmp = sdsnew (configs->out_dir_name);
    sds c = sdsfromlonglong (count);
    tmp = sdscat (tmp, "/V_t_");
    tmp = sdscat (tmp, c);
    FILE *f1 = fopen (tmp, "w");
    activity = print_grid_and_check_for_activity (the_grid, f1, count, save_in_binary);
    fclose (f1);
    sdsfree (tmp);
    sdsfree (c);
    return activity;
}

void set_spatial_stim(struct stim_config_hash *stim_configs, struct grid *the_grid) {

    struct stim_config *tmp = NULL;

    for (int i = 0; i < stim_configs->size; i++) {
        for (struct stim_config_elt *e = stim_configs->table[i % stim_configs->size]; e != 0; e = e->next) {
            tmp = e->value;
            tmp->set_spatial_stim (tmp, the_grid);
        }
    }
}

void set_ode_extra_data (struct extra_data_config *config, struct grid *the_grid, struct ode_solver *the_ode_solver) {

    free (the_ode_solver->edo_extra_data);
    the_ode_solver->edo_extra_data =
        config->set_extra_data (the_grid, config->config_data.config, &(the_ode_solver->extra_data_size));
}

void update_ode_state_vector (struct ode_solver *the_ode_solver, struct grid *the_grid, uint32_t max_number_of_cells) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    int n_edos = the_ode_solver->model_data.number_of_ode_equations;

    real *sv = the_ode_solver->sv;

	int i;

    if (the_ode_solver->gpu) {
#ifdef COMPILE_CUDA
        real *vms;
        size_t mem_size = max_number_of_cells * sizeof (real);

        vms = (real *)malloc (mem_size);
        check_cuda_errors (cudaMemcpy (vms, sv, mem_size, cudaMemcpyDeviceToHost));		

		#pragma omp parallel for
        for (i = 0; i < n_active; i++) {
            vms[ac[i]->sv_position] = (real)ac[i]->v;
        }

        check_cuda_errors (cudaMemcpy (sv, vms, mem_size, cudaMemcpyHostToDevice));
        free (vms);
#endif
    } else {
		#pragma omp parallel for
        for (i = 0; i < n_active; i++) {
            sv[ac[i]->sv_position * n_edos] = (real)ac[i]->v;
        }
    }
}

void save_old_cell_positions (struct grid *the_grid) 
{

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

	int i;

	#pragma omp parallel for
    for (i = 0; i < n_active; i++) 
    {
        ac[i]->sv_position = ac[i]->grid_position;
    }
}

void update_cells_to_solve (struct grid *the_grid, struct ode_solver *solver) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    if (solver->cells_to_solve) {
        free (solver->cells_to_solve);
    }

    solver->cells_to_solve = (uint32_t *)malloc (the_grid->num_active_cells * sizeof (uint32_t));
    uint32_t *cts = solver->cells_to_solve;
	int i;
	
	#pragma omp parallel for
    for (i = 0; i < n_active; i++) {
        cts[i] = ac[i]->sv_position;
    }
}

void set_initial_conditions_all_volumes (struct monodomain_solver *the_solver, struct grid *the_grid, double initial_v) 
{

    double alpha, h;
    struct cell_node **ac = the_grid->active_cells;
    uint32_t active_cells = the_grid->num_active_cells;
    double beta = the_solver->beta;
    double cm = the_solver->cm;
    double dt = the_solver->dt;
	int i;


	#pragma omp parallel for private(alpha, h)
    for (i = 0; i < active_cells; i++) 
    {
        h = ac[i]->face_length;
        alpha = ALPHA (beta, cm, dt, h);
        ac[i]->v = initial_v;
        ac[i]->b = initial_v * alpha;
    }
}

void update_monodomain (uint32_t initial_number_of_cells, uint32_t num_active_cells, struct cell_node **active_cells,
                        double beta, double cm, double dt_edp, real *sv, int n_equations_cell_model, bool use_gpu) 
{

    double h, alpha;

#ifdef COMPILE_CUDA
    real *vms = NULL;
    size_t mem_size = initial_number_of_cells * sizeof (real);

    if (use_gpu) 
    {
        vms = (real *)malloc (mem_size);
        check_cuda_errors (cudaMemcpy (vms, sv, mem_size, cudaMemcpyDeviceToHost));
    }
#endif
	int i;
	#pragma omp parallel for private(h, alpha)
    for (i = 0; i < num_active_cells; i++) 
    {
        h = active_cells[i]->face_length;
        alpha = ALPHA (beta, cm, dt_edp, h);

        if (use_gpu) 
        {
        #ifdef COMPILE_CUDA
            active_cells[i]->b = vms[active_cells[i]->sv_position] * alpha;
        #endif
        } 
        else 
        {
            active_cells[i]->b = sv[active_cells[i]->sv_position * n_equations_cell_model] * alpha;
        }
    }
#ifdef COMPILE_CUDA
    free (vms);
#endif
}

void print_solver_info (struct monodomain_solver *the_monodomain_solver, struct ode_solver *the_ode_solver,
                        struct grid *the_grid, struct user_options *options) {
    print_to_stdout_and_file ("System parameters: \n");
#if defined(_OPENMP)
    print_to_stdout_and_file ("Using OpenMP with %d threads\n", omp_get_max_threads ());
#endif
    if (the_ode_solver->gpu) {
        print_to_stdout_and_file ("Using GPU to solve ODEs\n");
    }

    print_to_stdout_and_file ("Initial V: %lf\n", the_ode_solver->model_data.initial_v);
    print_to_stdout_and_file ("Number of ODEs in cell model: %d\n", the_ode_solver->model_data.number_of_ode_equations);

    print_to_stdout_and_file ("Sigma X = %.10lf, Sigma Y = %.10lf, Sigma Z = %.10lf\n", the_monodomain_solver->sigma_x,
                              the_monodomain_solver->sigma_y, the_monodomain_solver->sigma_z);

    print_to_stdout_and_file ("Beta = %.10lf, Cm = %.10lf\n", the_monodomain_solver->beta, the_monodomain_solver->cm);

    print_to_stdout_and_file ("Initial N. of Elements = "
                              "%" PRIu32 "\n",
                              the_grid->num_active_cells);
    print_to_stdout_and_file ("PDE time step = %lf\n", the_monodomain_solver->dt);
    print_to_stdout_and_file ("ODE min time step = %lf\n", the_ode_solver->min_dt);
    print_to_stdout_and_file ("Simulation Final Time = %lf\n", the_monodomain_solver->final_time);
    print_to_stdout_and_file ("Maximum CG iterations = %d\n", the_monodomain_solver->max_iterations);
    print_to_stdout_and_file ("CG tolerance = %e\n", the_monodomain_solver->tolerance);
    if (the_monodomain_solver->use_jacobi) {
        print_to_stdout_and_file ("Using Jacobi preconditioner\n");
    }
    if (the_grid->adaptive) {
        print_to_stdout_and_file ("Using adaptativity\n");
        print_to_stdout_and_file ("Refinement Bound = %lf\n", the_monodomain_solver->refinement_bound);
        print_to_stdout_and_file ("Derefinement Bound = %lf\n", the_monodomain_solver->derefinement_bound);
        print_to_stdout_and_file ("Refining each %d time steps\n", the_monodomain_solver->refine_each);
        print_to_stdout_and_file ("Derefining each %d time steps\n", the_monodomain_solver->derefine_each);
    }

    print_to_stdout_and_file ("Print Rate = %d\n", options->print_rate);

    if (options->out_dir_name != NULL) {
        if (options->binary) {
            print_to_stdout_and_file ("Saving using binary output in %s dir\n", options->out_dir_name);

        } else {
            print_to_stdout_and_file ("Saving to plain text output in %s dir\n", options->out_dir_name);
        }
    } else {
        print_to_stdout_and_file ("The solution will not be saved\n");
    }

    if (options->stim_configs) 
    {
        print_to_stdout_and_file (LOG_LINE_SEPARATOR);

        if (options->stim_configs->size == 1)
            print_to_stdout_and_file ("Stimulus configuration:\n");
        else {
            print_to_stdout_and_file ("Stimuli configuration:\n");
        }

        for (int i = 0; i < options->stim_configs->size; i++) {
            for (struct stim_config_elt *e = options->stim_configs->table[i % options->stim_configs->size]; e != 0;
                 e = e->next) {

                print_to_stdout_and_file ("Stimulus name: %s\n", e->key);
                print_to_stdout_and_file ("Stimulus start: %lf\n", e->value->stim_start);
                print_to_stdout_and_file ("Stimulus duration: %lf\n", e->value->stim_duration);
                print_to_stdout_and_file ("Stimulus current: %lf\n", e->value->stim_current);
                print_to_stdout_and_file ("Stimulus library: %s\n", e->value->config_data.library_file_path);
                print_to_stdout_and_file ("Stimulus function: %s\n", e->value->config_data.function_name);
                struct string_hash *tmp = e->value->config_data.config;
                if (tmp->n == 1) {
                    print_to_stdout_and_file ("Stimulus extra parameter:\n");
                } else if (tmp->n > 1) {
                    print_to_stdout_and_file ("Stimulus extra parameters:\n");
                }

                STRING_HASH_PRINT_KEY_VALUE_LOG (tmp);

                print_to_stdout_and_file (LOG_LINE_SEPARATOR);
            }
        }
    }

    if (options->domain_config)
    {
        print_to_stdout_and_file ("Domain configuration:\n");
        print_to_stdout_and_file ("Domain name: %s\n", options->domain_config->domain_name);
        print_to_stdout_and_file ("Domain initial Space Discretization: %lf um\n", options->domain_config->start_h);

        if (the_grid->adaptive) {
            print_to_stdout_and_file ("Domain maximum Space Discretization: %lf um\n", options->domain_config->max_h);
            print_to_stdout_and_file ("The adaptivity will start in time: %lf ms\n",
                                    the_monodomain_solver->start_adapting_at);
        }

        if (options->domain_config->config_data.config->n == 1) {
            print_to_stdout_and_file ("Domain extra parameter:\n");
        } else if (options->domain_config->config_data.config->n > 1) {
            print_to_stdout_and_file ("Domain extra parameters:\n");
        }

        STRING_HASH_PRINT_KEY_VALUE_LOG (options->domain_config->config_data.config);
        print_to_stdout_and_file (LOG_LINE_SEPARATOR);
    }
    
    if (options->purkinje_config)
    {
        print_to_stdout_and_file ("Purkinje configuration:\n");
        print_to_stdout_and_file ("Purkinje network name: %s\n", options->purkinje_config->domain_name);
        print_to_stdout_and_file ("Purkinje network initial Space Discretization: %lf um\n", options->purkinje_config->start_h);

        if (options->purkinje_config->config_data.config->n == 1) 
        {
            print_to_stdout_and_file ("Purkinje extra parameter:\n");
        } 
        else if (options->purkinje_config->config_data.config->n > 1) {
            print_to_stdout_and_file ("Purkinje extra parameters:\n");
        }

        STRING_HASH_PRINT_KEY_VALUE_LOG (options->purkinje_config->config_data.config);
        print_to_stdout_and_file (LOG_LINE_SEPARATOR);
    }

    if (options->extra_data_config) {
        print_to_stdout_and_file ("Extra data ODE function configuration:\n");

        print_to_stdout_and_file ("Extra data library: %s\n",
                                  options->extra_data_config->config_data.library_file_path);
        print_to_stdout_and_file ("Extra data function: %s\n", options->extra_data_config->config_data.function_name);

        if (options->domain_config->config_data.config->n == 1) {
            print_to_stdout_and_file ("Extra data parameter:\n");
        } else if (options->domain_config->config_data.config->n > 1) {
            print_to_stdout_and_file ("Extra data parameters:\n");
        }

        STRING_HASH_PRINT_KEY_VALUE_LOG (options->extra_data_config->config_data.config);
        print_to_stdout_and_file (LOG_LINE_SEPARATOR);
    }
}

*/