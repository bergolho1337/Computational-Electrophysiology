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

#include "../grid/grid.h"
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
    
    the_monodomain_solver->M = nearbyint(options->final_time/options->dt);
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
    set_control_volumes(the_monodomain_solver,the_grid);
    set_derivative(the_monodomain_solver,the_grid);
    set_velocity_points(the_monodomain_solver,the_grid);
    set_plot_points(the_monodomain_solver);

    // Get a reference of the stimulus
    struct stim_config_hash *stimuli_configs = configs->stim_configs;

    double last_stimulus_time = -1.0;

    if (stimuli_configs) 
    {
        // Init all stimuli
        STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE_AND_KEY (stimuli_configs, init_stim_functions);

        //Find last stimuli
        size_t s_size = stimuli_configs->size;
        double s_end;
        for (int i = 0; i < (int)s_size; i++) 
        {
            for (struct stim_config_elt *e = stimuli_configs->table[i % s_size]; e != 0; e = e->next) 
            {
                s_end = e->value->stim_start + e->value->stim_duration;
                if(s_end > last_stimulus_time) last_stimulus_time = s_end;
            }
        }
    }

    if (the_monodomain_solver->use_steady_state)
        set_initial_conditions_from_file(the_monodomain_solver,configs);
    else
        set_initial_conditions_default(the_monodomain_solver);
    
    // Setting spatial stimulus
    if (stimuli_configs)
        set_spatial_stim(stimuli_configs, the_grid);

    // [DEBUG] Print out information about the simulation to be solved
    //print_solver_info (the_monodomain_solver, the_grid, configs);

    // Open the SST that will store the steady_state variables of the system on the time SST_RATE
    FILE *sst_file = fopen("output.sst","w+");

    int np = the_grid->the_purkinje_network->get_total_nodes();
    
    // Build the discretization matrix
    Eigen::SparseMatrix<double> A(np,np);
    set_matrix(A,the_monodomain_solver,the_grid);

    // Apply a LU decomposition over the matrix
    Eigen::SparseLU< Eigen::SparseMatrix<double> > sparseSolver(A);

    // Declare RHS and the solution vector
    Eigen::VectorXd b(np);
    Eigen::VectorXd x(np);

    std::cout << "[Solver] Time loop" << std::endl;
    int M = the_monodomain_solver->M;
    int print_rate = configs->print_rate;
    int sst_rate = configs->sst_rate;
    double dt = the_monodomain_solver->dt;
    

    // Time loop
    for (int i = 0; i < M; i++)
    {

        double t = i*dt;

        // Print the progress of the solution
        print_progress(i,M);
        
        // Write the solution data to a file
        write_plot_data(the_monodomain_solver,t);
        
        //if (i % 10 == 0)
        //    writeStateVectorToFile(i,max_num_digits,t);

        // Write the solution to .vtk file
        if (i % print_rate == 0) write_VTK_to_file(the_monodomain_solver,the_grid,i);
        if (i == sst_rate-1) write_steady_state_to_file(sst_file,the_monodomain_solver);
        
        // Solve the PDE (diffusion phase)
        assemble_load_vector(b,the_monodomain_solver,the_grid);
        x = sparseSolver.solve(b);

        move_v_star(x,the_monodomain_solver);

        // Solve the ODEs (reaction phase)
        solve_odes(t,the_monodomain_solver,configs->stim_configs);

        // Calculate the maximum derivative for each volume
        //calc_max_derivative(t,stim_config->start_period);

        // Jump to the next iteration
        next_timestep(the_monodomain_solver);
    }

}

double* merge_stimulus (struct stim_config_hash *stim_configs,\
                    const int np, const double cur_time)
{
    double *merged_stims = (double*)calloc(sizeof(double),np);
    
    struct stim_config *tmp = NULL;
    double time = cur_time;

    double stim_period;
    double stim_start, stim_duration;
    double start_period, end_period, period_step;
    int n_cycles;

    double new_time;

	int i;
    if(stim_configs) 
    {
        for (int k = 0; k < (int)stim_configs->size; k++) 
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
                    //#pragma omp parallel for
                    for (i = 0; i < np; i++) 
                    {
                        merged_stims[i] = tmp->spatial_stim_currents[i];
                    }
                }
                time = cur_time;
            }
        }
    }

    return merged_stims;
}

void solve_odes (const double t,\
                 struct monodomain_solver *solver,
                 struct stim_config_hash *stim_configs)
{
    double *merged_stims = merge_stimulus(stim_configs,solver->num_volumes,t);

    solve_model_ode_cpu_fn *solve_odes_pt = solver->solve_model_ode_cpu;
    solve_odes_pt(solver->dt,merged_stims,solver->num_volumes,solver->volumes);

    free(merged_stims);
}

void set_spatial_stim(struct stim_config_hash *stim_configs,\
                      struct grid *the_grid) 
{
    struct stim_config *tmp = NULL;

    for (int i = 0; i < stim_configs->size; i++) 
    {
        for (struct stim_config_elt *e = stim_configs->table[i % stim_configs->size]; e != 0; e = e->next) 
        {
            tmp = e->value;
            tmp->set_spatial_stim (tmp, the_grid);
        }
    }
}

void set_initial_conditions_from_file (struct monodomain_solver *solver, struct user_options *options)
{
    std::cout << "[Solver] Setting initial conditions from SST file: " << options->sst_filename << std::endl;
    
    // Initialize the solver with the 'initial_v' and 'number_quations' from the celular model
    (*(solver->get_cell_model_data))(&(solver->model_data));
    int n_odes = solver->model_data.number_of_ode_equations;

    FILE *sst_file = fopen(options->sst_filename,"r");
    if (!sst_file)
    {
        cerr << "[-] ERROR! Could not open SST file!" << endl;
        exit(EXIT_FAILURE);
    }

    int np = solver->num_volumes;
    // Iterate over all the control volumes
    for (int i = 0; i < np; i++)
    {
        // Iterate over all the ODEs equations
        for (int j = 0; j < n_odes; j++)
        {
            if (!fscanf(sst_file,"%lf",&solver->volumes[i].y_old[j]))
            {
                cerr << "[-] ERROR! Could not open SST file!" << endl;
                exit(EXIT_FAILURE); 
            }
        }
    }
    fclose(sst_file);
}

void set_initial_conditions_default (struct monodomain_solver *solver)
{
    std::cout << "[Solver] Setting default initial conditions" << std::endl;

    // Initialize the solver with the 'initial_v' and 'number_quations' from the celular model
    (*(solver->get_cell_model_data))(&(solver->model_data));
    int n_odes = solver->model_data.number_of_ode_equations;

    set_ode_initial_conditions_cpu_fn *soicc_fn_pt = solver->set_ode_initial_conditions_cpu;

    if(!soicc_fn_pt) 
    {
        fprintf(stderr, "The ode solver was set to use the CPU, \n "
                "but no function called set_model_initial_conditions_cpu "
                "was provided in the %s shared library file\n", solver->model_data.model_library_path);
        exit(11);
    }
    //#pragma omp parallel for
    for(int i = 0; i < solver->num_volumes; i++) 
    {
        soicc_fn_pt(solver->volumes[i].y_old);
    }
}

void set_matrix (Eigen::SparseMatrix<double> &a,\
                struct monodomain_solver *solver,\
                struct grid *the_grid)
{
    std::cout << "[Solver] Building discretization matrix" << std::endl;

    double g_GAP = solver->G_gap;
    double diameter = solver->start_diameter;
    double sigma_c = solver->sigma_c;
    double beta = solver->beta;
    double cm = solver->cm;
    double dt = solver->dt;
    double dx = the_grid->dx;

    // Compute the coefficients values
    double A = (4.0*g_GAP*dx) / (M_PI*diameter*diameter);
    double B = (sigma_c);
    double C = (beta*cm*dx*dx) / (dt);

    // Non-zero coefficients
    vector< Eigen::Triplet<double> > coeff;

    double diagonal_value;
    Node *ptr = the_grid->the_purkinje_network->get_list_nodes();
    while (ptr != NULL)
    {
        int u = ptr->id;
        Edge *ptrl = ptr->list_edges;
        diagonal_value = C;

        while (ptrl != NULL)
        {
            double value;
            int v = ptrl->id;
            int link_type = ptrl->link_type;

            // Citoplasm link
            if (link_type == 0)
            {
                value = -B;
                diagonal_value += B;
            }
            // Gap junction link
            else
            {
                value = -A;
                diagonal_value += A;
            }
            coeff.push_back(Eigen::Triplet<double>(u,v,value));

            ptrl = ptrl->next;
        }
        coeff.push_back(Eigen::Triplet<double>(u,u,diagonal_value));

        ptr = ptr->next;
    }

    a.setFromTriplets(coeff.begin(),coeff.end());
    a.makeCompressed();
}

void assemble_load_vector (Eigen::VectorXd &b,\
                            struct monodomain_solver *solver,\
                            struct grid *the_grid)
{
    double beta = solver->beta;
    double cm = solver->cm;
    double dt = solver->dt;
    double dx = the_grid->dx;

    double C = (beta*cm*dx*dx) / (dt);

    int np = b.size();
    //#pragma omp parallel for
    for (int i = 0; i < np; i++)
        b(i) = solver->volumes[i].y_old[0] * C;

}

void move_v_star (const Eigen::VectorXd vm,\
                    struct monodomain_solver *solver)
{
    int neq = solver->model_data.number_of_ode_equations;
    int np = vm.size();
    //#pragma omp parallel for
    for (int i = 0; i < np; i++)
    {
        solver->volumes[i].y_star[0] = vm(i);
        for (int j = 1; j < neq; j++)
            solver->volumes[i].y_star[j] = solver->volumes[i].y_old[j];
    }
}

void swap (double **a, double **b)
{
    double *tmp = *a;
    *a = *b;
    *b = tmp;
}

void next_timestep (struct monodomain_solver *solver)
{
    int np = solver->num_volumes;
    //#pragma omp parallel for
    for (int i = 0; i < np; i++) 
        swap(&solver->volumes[i].y_old,&solver->volumes[i].y_new);
}

void print_solver_info (struct monodomain_solver *the_monodomain_solver,\
                        struct grid *the_grid,\
                        struct user_options *configs)
{
    std::cout << "////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "number of threads = " << the_monodomain_solver->num_threads << std::endl;
    std::cout << "dt = " << the_monodomain_solver->dt << std::endl;
    std::cout << "tmax = " << the_monodomain_solver->final_time << std::endl;
    std::cout << "number of timesteps = " << the_monodomain_solver->M << std::endl;
    std::cout << "dx = " << the_grid->dx << std::endl;
    std::cout << "network_filename = " << configs->network_filename << std::endl;
    std::cout << "steady_state_filename = " << configs->sst_filename << std::endl;
    std::cout << "plot_filename = " << configs->plot_filename << std::endl;
    std::cout << "print_rate = " << configs->print_rate << std::endl;
    std::cout << "sst_rate = " << configs->sst_rate << std::endl;
    std::cout << "diameter = " << configs->start_diameter << std::endl;
    std::cout << "beta = " << the_monodomain_solver->beta << std::endl;
    std::cout << "Cm = " << the_monodomain_solver->cm << std::endl;
    std::cout << "sigma_c = " << the_monodomain_solver->sigma_c << std::endl;
    std::cout << "G_gap = " << the_monodomain_solver->G_gap << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////" << std::endl;
    
    if (configs->stim_configs) 
    {

        if (configs->stim_configs->size == 1)
            std::cout << "Stimulus configuration:" << std::endl;
        else
            std::cout << "Stimuli configuration:" << std::endl;

        for (int i = 0; i < (int)configs->stim_configs->size; i++) 
        {
            for (struct stim_config_elt *e = configs->stim_configs->table[i % configs->stim_configs->size]; e != 0;
                 e = e->next) 
            {

                std::cout << "----------------------------------------------------------------------------" << std::endl;
                std::cout << "Stimulus name: " << e->key << std::endl;
                std::cout << "Stimulus start: " << e->value->stim_start << std::endl;
                std::cout << "Stimulus duration: " << e->value->stim_duration << std::endl;
                std::cout << "Stimulus current: " << e->value->stim_current << std::endl;
                std::cout << "Stimulus start period: " << e->value->start_period << std::endl;
                std::cout << "Stimulus end period: " << e->value->end_period << std::endl;
                std::cout << "Stimulus period step: " << e->value->period_step << std::endl;
                std::cout << "Stimulus function: " << e->value->config_data.function_name << std::endl;
                std::cout << "Number of cycles: " << e->value->n_cycles << std::endl;
                struct string_hash *tmp = e->value->config_data.config;
                if (tmp->n == 1) 
                    std::cout << "Stimulus extra parameter:" << std::endl;
                else if (tmp->n > 1) 
                    std::cout << "Stimulus extra parameters:" << std::endl;

                STRING_HASH_PRINT_KEY_VALUE_LOG (tmp);

            }
        }
    }
}

void print_progress (int iter, int max_iter)
{
    double progress = iter / (double)max_iter;
    
    cout << "Progress: " << int(progress * 100.0) << " %\r";
    cout.flush();
}

void write_VTK_to_file (struct monodomain_solver *solver,\
                        struct grid *the_grid, int iter)
{
    FILE *file;
    int np, ne;
    char filename[50];
    Node *ptr = the_grid->the_purkinje_network->get_list_nodes();
    np = the_grid->the_purkinje_network->get_total_nodes();
    ne = the_grid->the_purkinje_network->get_total_edges();

    // Write the transmembrane potential
    sprintf(filename,"vtk/sol%d.vtk",iter);
    file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Monodomain MVF\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %d float\n",np);
    while (ptr != NULL)
    {
        fprintf(file,"%e %e %e\n",ptr->x,ptr->y,ptr->z);
        ptr = ptr->next;
    }
    fprintf(file,"LINES %d %d\n",ne,ne*3);
    ptr = the_grid->the_purkinje_network->get_list_nodes();
    while (ptr != NULL)
    {
        Edge *ptrl = ptr->list_edges;
        while (ptrl != NULL)
        {
            fprintf(file,"2 %d %d\n",ptr->id,ptrl->dest->id);
            ptrl = ptrl->next;
        }
        ptr = ptr->next;
    }

    fprintf(file,"POINT_DATA %d\n",np);
    fprintf(file,"SCALARS vm float 1\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    ptr = the_grid->the_purkinje_network->get_list_nodes();
    while (ptr != NULL)
    {
        fprintf(file,"%e\n",solver->volumes[ptr->id].y_old[0]);
        ptr = ptr->next;
    }
    fclose(file);
}

void write_plot_data (struct monodomain_solver *solver, double t)
{
    for (int i = 1; i < solver->plot->np; i++)
        fprintf(solver->plot->plotFile[i-1],"%.10lf %.10lf\n",t,solver->volumes[solver->plot->ids[i]].y_old[0]);
}

void write_steady_state_to_file (FILE *sst_file,\
                                struct monodomain_solver *solver)
{
    int neq = solver->model_data.number_of_ode_equations;
    int np = solver->num_volumes;
    for (int i = 0; i < np; i++)
    {
        for (int j = 0; j < neq; j++)
            fprintf(sst_file,"%.10lf ",solver->volumes[i].y_old[j]);
        fprintf(sst_file,"\n");
    }
}
