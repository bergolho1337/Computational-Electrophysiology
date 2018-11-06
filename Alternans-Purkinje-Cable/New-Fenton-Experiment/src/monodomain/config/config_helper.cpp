#include "config_helper.h"

void set_celular_model (struct monodomain_solver *solver, struct user_options *configs)
{
    void *handle = dlopen (configs->model_file_path, RTLD_LAZY);
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
    std::cout << "Library " << configs->model_file_path <<  " loaded with sucess" << std::endl;
}

void set_control_volumes (struct monodomain_solver *solver, struct grid *the_grid)
{
    // Copy the number of control volumes to the 'monodomain_solver' structure
    solver->num_volumes = the_grid->the_purkinje_network->get_total_nodes();

    // Capture the number of equations of the celullar model
    int neq = solver->model_data.number_of_ode_equations;
    int np = the_grid->the_purkinje_network->get_total_nodes();
    
    // Allocate memory for the structure
    solver->volumes = (struct control_volume*)malloc(sizeof(struct control_volume)*np);
    for (int i = 0; i < np; i++)
    {
        solver->volumes[i].y_old = (double*)calloc(neq,sizeof(double));
        solver->volumes[i].y_new = (double*)calloc(neq,sizeof(double));
        solver->volumes[i].y_star = (double*)calloc(neq,sizeof(double));
    }
}

void set_derivative (struct monodomain_solver *solver, struct grid *the_grid)
{
    int np = the_grid->the_purkinje_network->get_total_nodes();
    solver->dvdt = (struct derivative*)malloc(np*sizeof(struct derivative));
    for (int i = 0; i < np; i++) 
        solver->dvdt[i].value = 0;
}

void set_velocity_points (struct monodomain_solver *solver, struct grid *the_grid)
{
    solver->vel = (struct velocity*)malloc(sizeof(struct velocity));
    solver->vel->velocity_file = fopen("output/velocity.txt","w+");

    // First point is always the source
    solver->vel->np = solver->plot->np-1;
    solver->vel->id_source = solver->plot->ids[0];
    solver->vel->ids = (int*)malloc(sizeof(int)*solver->vel->np);
    solver->vel->t2 = (double*)malloc(sizeof(double)*solver->vel->np);
    for (int i = 0; i < solver->vel->np; i++) 
        solver->vel->ids[i] = solver->plot->ids[i+1];
}

void set_plot_points (struct monodomain_solver *solver)
{
    char filename[200];
    solver->plot->plotFile = (FILE**)malloc(sizeof(FILE*)*(solver->plot->np-1));
    for (int i = 1; i < solver->plot->np; i++)
    {
        solver->plot->plotFile[i-1] = (FILE*)malloc(sizeof(FILE));
        sprintf(filename,"output/data-%d.dat",solver->plot->ids[i]);
        solver->plot->plotFile[i-1] = fopen(filename,"w+");
    }
}