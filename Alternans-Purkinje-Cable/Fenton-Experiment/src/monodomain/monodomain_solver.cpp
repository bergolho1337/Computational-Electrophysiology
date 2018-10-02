#include "monodomain_solver.h"
#include "ode_solver.h"
#include <cassert>
#include <cinttypes>
#include <vector>
#include <string>
#include <sstream>
#include "../utils/logfile_utils.h"
#include "../utils/stop_watch.h"
#include "../purkinje/purkinje.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

//#include "config/purkinje_config.h"

static inline double ALPHA (double beta, double cm, double dt, double h) 
{
    return (((beta * cm) / dt) * UM2_TO_CM2) * pow (h, 3.0);
    //return (((beta * cm) / dt)) * pow (h, 3.0);
}

struct monodomain_solver *new_monodomain_solver ()
{
    struct monodomain_solver *result = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));

    result->beta = 0.14;
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

    the_monodomain_solver->dt = options->dt;

    the_monodomain_solver->sigma = options->sigma;
    the_monodomain_solver->beta = options->beta;
    the_monodomain_solver->cm = options->cm;
}

void solve_monodomain (struct monodomain_solver *monodomain_solver, struct ode_solver *ode_solver,\
                        struct grid *grid, struct user_options *configs)
{
    assert (configs);
    assert (grid);
    assert (monodomain_solver);
    assert (ode_solver);

    print_to_stdout_and_file (LOG_LINE_SEPARATOR);

    long ode_total_time = 0, lu_total_time = 0, total_write_time = 0, total_mat_time = 0, total_ref_time = 0,
         total_deref_time = 0, lu_partial, total_config_time = 0;

    struct stop_watch solver_time, ode_time, lu_time, part_solver, part_mat, write_time, ref_time, deref_time,
        config_time;

    init_stop_watch (&config_time);

    start_stop_watch (&config_time);

    ///////MAIN CONFIGURATION BEGIN//////////////////
    init_ode_solver_with_cell_model(ode_solver);
    struct stim_config_hash *stimuli_configs = configs->stim_configs;
    struct purkinje_config *pk_config = configs->purkinje_config;

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
    if (pk_config) 
    {
        set_spatial_purkinje(pk_config,grid);
    } 
    else 
    {
        print_to_stdout_and_file ("No Purkinje configuration provided! Exiting!\n");
        exit (EXIT_FAILURE);
    }

    bool gpu = ode_solver->gpu;
    int count = 0;

    bool save_to_file = (configs->out_dir_name != NULL);

    double dt_edp = monodomain_solver->dt;
    double finalT = monodomain_solver->final_time;

    double beta = monodomain_solver->beta;
    double cm = monodomain_solver->cm;

    double dt_edo = ode_solver->dt;

    #ifdef COMPILE_CUDA
    if (gpu) 
    {
        int device_count;
        int device = ode_solver->gpu_id;
        check_cuda_errors (cudaGetDeviceCount (&device_count));
        struct cudaDeviceProp prop;
        check_cuda_errors (cudaGetDeviceProperties (&prop, ode_solver->gpu_id));
        print_to_stdout_and_file ("%d devices available, running on Device %d: %s\n", device_count, device, prop.name);
        check_cuda_errors (cudaSetDevice (device));
    }
    #endif
    
    order_grid_cells (grid);
    uint32_t original_num_cells = grid->num_active_cells;

    save_old_cell_positions (grid);

    print_to_stdout_and_file ("Setting ODE's initial conditions\n");
    set_ode_initial_conditions_for_all_volumes (ode_solver, grid->num_active_cells);

    double initial_v = ode_solver->model_data.initial_v;

    total_config_time = stop_stop_watch (&config_time);

    print_solver_info (monodomain_solver, ode_solver, grid, configs);

    int ode_step = 1;

    fflush (stdout);

    init_stop_watch (&solver_time);
    init_stop_watch (&ode_time);
    init_stop_watch (&lu_time);
    init_stop_watch (&part_solver);
    init_stop_watch (&part_mat);
    init_stop_watch (&write_time);
    init_stop_watch (&ref_time);
    init_stop_watch (&deref_time);

    print_to_stdout_and_file ("Assembling Monodomain Matrix Begin\n");
    start_stop_watch (&part_mat);

    set_initial_conditions_all_volumes (monodomain_solver, grid, initial_v);

    // Assemble the matrix of the PDE
    Eigen::SparseMatrix<double> A = assembly_matrix(monodomain_solver,grid,pk_config);
    // Sparse LU Decomposition
    Eigen::SparseLU< Eigen::SparseMatrix<double> > sparse_solver(A);

    total_mat_time = stop_stop_watch (&part_mat);
    print_to_stdout_and_file ("Assembling Monodomain Matrix End\n");
    print_to_stdout_and_file (LOG_LINE_SEPARATOR);

    start_stop_watch (&solver_time);

    int print_rate = configs->print_rate;

    double solver_error;
    uint32_t solver_iterations = 0;

    if (stimuli_configs)
        set_spatial_stim(stimuli_configs, grid);

    double cur_time = 0.0;
    bool activity;
    bool save_in_binary = false;

    // Declare RHS and the solution vector
    Eigen::VectorXd b(original_num_cells);
    Eigen::VectorXd x(original_num_cells);

    print_to_stdout_and_file ("Starting simulation\n");

    while (cur_time <= finalT)
    {
        if (save_to_file)
        {
            if (count % print_rate == 0)
            {
                start_stop_watch (&write_time);

                activity = print_result(grid, configs, count, save_in_binary);

                total_write_time += stop_stop_watch (&write_time);
            }

            if (count % 60000 == 0)
            {
                //print_steady_state(grid,ode_solver,configs,count);
            }
        }

        if (cur_time > 0.0) 
        {
            // Pass the solution to next iteration: n+1 = n
            update_ode_state_vector (ode_solver, grid, original_num_cells);
        }

        start_stop_watch (&ode_time);

        // Solve the ODE for all volumes
        solve_all_volumes_odes (ode_solver, grid->num_active_cells, cur_time, ode_step, stimuli_configs);

        // Update the 'b' value of each cell (Remember, the 'b' variable is the RHS of the linear system)
        update_monodomain (original_num_cells, grid->num_active_cells, grid->active_cells, beta, cm, dt_edp,
                           ode_solver->sv, ode_solver->model_data.number_of_ode_equations, gpu);

        ode_total_time += stop_stop_watch (&ode_time);

        // Solve the linear system of the PDE
        start_stop_watch (&lu_time);

        assembly_load_vector(b,grid->active_cells,grid->num_active_cells);
        // Sparse Back+Forward Substitution
        x = sparse_solver.solve(b);
        move_solution_to_cells(x,grid->active_cells,grid->num_active_cells);

        lu_partial = stop_stop_watch (&lu_time);

        lu_total_time += lu_partial;

        if (count % print_rate == 0) 
        {
            print_to_stdout_and_file ("t = %lf, Number of Cells: %u\n",cur_time,grid->num_active_cells);
        }

        count++;
        cur_time += dt_edp;

    }
    print_to_stdout_and_file ("Resolution Time: %ld μs\n", stop_stop_watch (&solver_time));
    print_to_stdout_and_file ("ODE Total Time: %ld μs\n", ode_total_time);
    print_to_stdout_and_file ("LU Total Time: %ld μs\n", lu_total_time);
    print_to_stdout_and_file ("Mat time: %ld μs\n", total_mat_time);
    print_to_stdout_and_file ("Write time: %ld μs\n", total_write_time);
    print_to_stdout_and_file ("Initial configuration time: %ld μs\n", total_config_time);

}

Eigen::SparseMatrix<double> assembly_matrix (struct monodomain_solver *monodomain_solver, struct grid *grid, struct purkinje_config *pk_config) 
{
    print_to_stdout_and_file("Building Eigen Sparse LU Matrix\n");

    struct graph *pk_net = grid->the_purkinje_network;
    int n = pk_net->total_nodes;
    double cm = monodomain_solver->cm;
    double beta = monodomain_solver->beta;
    double sigma = monodomain_solver->sigma;
    double dt = monodomain_solver->dt;
    double h = pk_config->start_h;

    Eigen::SparseMatrix<double> A(n,n);
    std::vector< Eigen::Triplet<double> > coeff;

    double alpha = ALPHA(beta,cm,dt,h);
    
    struct node *ptr = pk_net->list_nodes;
    while (ptr != NULL)
    {
        int u = ptr->id;
        double value = -sigma * h;
        
        struct edge *ptrl = ptr->list_edges;
        while (ptrl != NULL)
        {
            int v = ptrl->dest->id;
            coeff.push_back( Eigen::Triplet<double>(u,v,value) );
            ptrl = ptrl->next;
        }
        value = (ptr->num_edges * sigma * h) + alpha;
        coeff.push_back( Eigen::Triplet<double>(u,u,value) );
        ptr = ptr->next;
    }

    // Print non-zero coefficients
    //FILE *file = fopen("matrix.txt","w+");
    //for (int i = 0; i < coeff.size(); i++)
    //    fprintf(file,"(%d,%d) = %.10lf\n",coeff[i].row(),coeff[i].col(),coeff[i].value());
    //fclose(file);

    A.setFromTriplets(coeff.begin(),coeff.end());
    A.makeCompressed();

    return A;
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

void set_spatial_stim(struct stim_config_hash *stim_configs, struct grid *the_grid) 
{
    struct stim_config *tmp = NULL;

    for (int i = 0; i < stim_configs->size; i++) 
    {
        for (struct stim_config_elt *e = stim_configs->table[i % stim_configs->size]; e != 0; e = e->next) 
        {
            tmp = e->value;
            set_stimulus(tmp, the_grid);
        }
    }
}

void set_stimulus (struct stim_config *config, struct grid *grid)
{
    uint32_t n_active = grid->num_active_cells;
    struct cell_node **ac = grid->active_cells;

    bool stim;
    float stim_current = config->stim_current;
    float stim_value;

    int id;
    id = config->id_limit;

    if(config->spatial_stim_currents) 
    {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (float*)malloc(n_active*sizeof(float));

	int i;

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) 
    {
        stim = i <= id;

        if (stim) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}

// Move the solution to the next step
void update_ode_state_vector (struct ode_solver *the_ode_solver, struct grid *the_grid, uint32_t max_number_of_cells) 
{

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    int n_edos = the_ode_solver->model_data.number_of_ode_equations;

    float *sv = the_ode_solver->sv;

	int i;

    if (the_ode_solver->gpu) 
    {
#ifdef COMPILE_CUDA
        float *vms;
        size_t mem_size = max_number_of_cells * sizeof (float);

        vms = (float *)malloc (mem_size);
        check_cuda_errors (cudaMemcpy (vms, sv, mem_size, cudaMemcpyDeviceToHost));		

		#pragma omp parallel for
        for (i = 0; i < n_active; i++) 
        {
            vms[ac[i]->sv_position] = (float)ac[i]->v;
        }

        check_cuda_errors (cudaMemcpy (sv, vms, mem_size, cudaMemcpyHostToDevice));
        free (vms);
#endif
    } 
    else 
    {
		#pragma omp parallel for
        for (i = 0; i < n_active; i++) 
        {
            sv[ac[i]->sv_position * n_edos] = (float)ac[i]->v;
        }
    }
}

void update_monodomain (uint32_t initial_number_of_cells, uint32_t num_active_cells, struct cell_node **active_cells,
                        double beta, double cm, double dt_edp, float *sv, int n_equations_cell_model, bool use_gpu) 
{

    double h, alpha;

#ifdef COMPILE_CUDA
    float *vms = NULL;
    size_t mem_size = initial_number_of_cells * sizeof (float);

    if (use_gpu) 
    {
        vms = (float *)malloc (mem_size);
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
                        struct grid *the_grid, struct user_options *options) 
                        {
    print_to_stdout_and_file ("System parameters: \n");
#if defined(_OPENMP)
    print_to_stdout_and_file ("Using OpenMP with %d threads\n", omp_get_max_threads ());
#endif
    if (the_ode_solver->gpu) {
        print_to_stdout_and_file ("Using GPU to solve ODEs\n");
    }

    print_to_stdout_and_file ("Initial V: %lf\n", the_ode_solver->model_data.initial_v);
    print_to_stdout_and_file ("Number of ODEs in cell model: %d\n", the_ode_solver->model_data.number_of_ode_equations);

    print_to_stdout_and_file ("Sigma X = %.10lf\n", the_monodomain_solver->sigma);

    print_to_stdout_and_file ("Beta = %.10lf, Cm = %.10lf\n", the_monodomain_solver->beta, the_monodomain_solver->cm);

    print_to_stdout_and_file ("Initial N. of Elements = "
                              "%" PRIu32 "\n",
                              the_grid->num_active_cells);
    print_to_stdout_and_file ("PDE time step = %lf\n", the_monodomain_solver->dt);
    print_to_stdout_and_file ("Simulation Final Time = %lf\n", the_monodomain_solver->final_time);

    print_to_stdout_and_file ("Print Rate = %d\n", options->print_rate);

    if (options->out_dir_name != NULL) {
        print_to_stdout_and_file ("Saving to plain text output in %s dir\n", options->out_dir_name);
    } else {
        print_to_stdout_and_file ("The solution will not be saved\n");
    }
    if (options->out_steady_state_dir != NULL) {
        print_to_stdout_and_file ("Saving steady-state solution to plain text output in %s dir\n", options->out_steady_state_dir);
    } else {
        print_to_stdout_and_file ("The solution will not be saved\n");
    }

    if (options->stim_configs) 
    {
        print_to_stdout_and_file (LOG_LINE_SEPARATOR);

        print_to_stdout_and_file ("Stimulus configuration:\n");

        for (int i = 0; i < options->stim_configs->size; i++) 
        {
            for (struct stim_config_elt *e = options->stim_configs->table[i % options->stim_configs->size]; e != 0;
                 e = e->next) 
                 {

                print_to_stdout_and_file ("Stimulus name: %s\n", e->key);
                print_to_stdout_and_file ("Stimulus start: %lf\n", e->value->stim_start);
                print_to_stdout_and_file ("Stimulus duration: %lf\n", e->value->stim_duration);
                print_to_stdout_and_file ("Stimulus current: %lf\n", e->value->stim_current);
                print_to_stdout_and_file ("Stimulus start period: %lf\n", e->value->start_period);
                print_to_stdout_and_file ("Stimulus end period: %lf\n", e->value->end_period);
                print_to_stdout_and_file ("Stimulus period step: %lf\n", e->value->period_step);
                print_to_stdout_and_file ("Stimulus number of cycles: %d\n", e->value->n_cycles);
                print_to_stdout_and_file ("Stimulus id limit: %d\n", e->value->id_limit);

                print_to_stdout_and_file (LOG_LINE_SEPARATOR);
            }
        }

        print_to_stdout_and_file (LOG_LINE_SEPARATOR);

    }
    
    if (options->purkinje_config)
    {
        print_to_stdout_and_file ("Purkinje configuration:\n");
        print_to_stdout_and_file ("Purkinje network name: %s\n", options->purkinje_config->network_name);
        print_to_stdout_and_file ("Purkinje network initial Space Discretization: %lf um\n", options->purkinje_config->start_h);
        print_to_stdout_and_file ("Purkinje network initial Diameter: %lf um\n", options->purkinje_config->diameter);
        print_to_stdout_and_file ("Purkinje network filename: %s\n", options->purkinje_config->network_filename);

        print_to_stdout_and_file (LOG_LINE_SEPARATOR);
    }

}

bool print_result(const struct grid *the_grid, const struct user_options *configs, int count, bool save_in_binary) 
{
    bool activity = true;
    std::stringstream ss;
    ss << configs->out_dir_name << "/V_t_" << count; 
    std::string tmp = ss.str();

    FILE *f1 = fopen (tmp.c_str(), "w");
    
    activity = print_grid_and_check_for_activity (the_grid, f1, count, save_in_binary);
    
    fclose (f1);
    
    return activity;
}

void assembly_load_vector (Eigen::VectorXd &b, struct cell_node **active_cells, uint32_t num_active_cells)
{
    struct cell_node **ac = active_cells;
    #pragma omp parallel for
    for (int i = 0; i < num_active_cells; i++)
    {
        b(i) = ac[i]->b;
    }
}

void move_solution_to_cells (Eigen::VectorXd &x, struct cell_node **active_cells, uint32_t num_active_cells)
{
    struct cell_node **ac = active_cells;
    #pragma omp parallel for
    for (int i = 0; i < num_active_cells; i++)
    {
        ac[i]->v = x(i);
    }
}
