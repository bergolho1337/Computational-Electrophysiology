#include "monodomain_solver.h"

struct monodomain_solver* new_monodomain_solver ()
{
    struct monodomain_solver *result = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));
    
    // Default values
    result->cm = 1.2;
    result->beta = 0.14;
    result->sigma_c = 0.004;
    result->G_gap = 0.628;
    
    return result;
}

void configure_monodomain_from_options (struct monodomain_solver *solver, struct user_options *options)
{
    assert(solver);
    assert(options);

    double dt = options->dt;
    solver->dt = dt;

    double tmax = options->tmax;
    solver->tmax = tmax;

    int use_steady_state = options->use_steady_state;
    if (use_steady_state)
        solver->use_steady_state = true;
    else
        solver->use_steady_state = false;

    // Calculate the beta parameter based on diameter of the Purkinje cell
    double diameter = options->diameter;
    solver->beta = 4.0f / diameter * 1.0e-04;

    double sigma_c = options->sigma_c;
    solver->sigma_c = sigma_c;

    double G_gap = options->G_gap;
    solver->G_gap = G_gap;

}

void print_solver_info (struct monodomain_solver *the_monodomain_solver,\
                        struct ode_solver *the_ode_solver,\
                        struct graph *the_purkinje_network,\
                        struct user_options *options)
{
    assert(the_monodomain_solver);
    assert(the_ode_solver);
    assert(the_purkinje_network);
    assert(options);

    fprintf(stdout,"%s\n",PRINT_LINE2);
    fprintf(stdout,"[Monodomain] Beta = %.10lf\n",the_monodomain_solver->beta);
    fprintf(stdout,"[Monodomain] Cm = %.10lf\n",the_monodomain_solver->cm);
    fprintf(stdout,"[Monodomain] sigma_c = %.10lf\n",the_monodomain_solver->sigma_c);
    fprintf(stdout,"[Monodomain] G_gap = %.10lf\n",the_monodomain_solver->G_gap);
    fprintf(stdout,"[Monodomain] dt = %.10lf\n",the_monodomain_solver->dt);
    fprintf(stdout,"[Monodomain] tmax = %.10lf\n",the_monodomain_solver->tmax);
    fprintf(stdout,"[Monodomain] use_steady_state = %d\n",the_monodomain_solver->use_steady_state);
    fprintf(stdout,"%s\n",PRINT_LINE);
    fprintf(stdout,"[ODE] Celular model name = %s\n",the_ode_solver->model_name);
    fprintf(stdout,"[ODE] Number of ODE equations = %d\n",the_ode_solver->num_ode_equations);
    fprintf(stdout,"[ODE] Number of volumes = %d\n",the_ode_solver->n_active_cells);
    fprintf(stdout,"%s\n",PRINT_LINE);
    fprintf(stdout,"[Purkinje] Total nodes = %d\n",the_purkinje_network->total_nodes);
    fprintf(stdout,"[Purkinje] Total edges = %d\n",the_purkinje_network->total_edges);
    fprintf(stdout,"[Purkinje] Diameter = %.10lf\n",the_purkinje_network->diameter);
    fprintf(stdout,"[Purkinje] Spatial discretization = %.10lf\n",the_purkinje_network->dx);
    fprintf(stdout,"[Purkinje] Number of cell divisions = %d\n",the_purkinje_network->num_div_cell);
    fprintf(stdout,"%s\n",PRINT_LINE);
    fprintf(stdout,"[Options] Network filename = %s\n",options->pk_network_filename);
    fprintf(stdout,"[Options] Plot points filename = %s\n",options->plot_filename);
    fprintf(stdout,"%s\n",PRINT_LINE2);

}

Eigen::SparseMatrix<double> assemble_matrix (struct monodomain_solver *the_monodomain_solver,\
                      struct graph *the_purkinje_network)
{
    assert(the_monodomain_solver);
    assert(the_purkinje_network);

    uint32_t nc = the_purkinje_network->total_nodes;

    double dt = the_monodomain_solver->dt; 
    double beta = the_monodomain_solver->beta;
    double cm = the_monodomain_solver->cm;
    double sigma_c = the_monodomain_solver->sigma_c;
    double G_gap = the_monodomain_solver->G_gap;
    double h = the_purkinje_network->dx;
    double d = the_purkinje_network->diameter;

    //double A = (beta*cm*h*h)/(dt);
    //double B = sigma_c;
    //double C = (4.0*G_gap)/(M_PI*d*d);

    double A = (beta*cm*h*h*h) / dt;
    double B = sigma_c*h;
    double C = sigma_c*h;

    Eigen::SparseMatrix<double> a(nc,nc);
    std::vector< Eigen::Triplet<double> > coeff;

    struct node *ptr = the_purkinje_network->list_nodes;
    while (ptr != NULL)
    {
        double value;
        uint32_t u = ptr->id;
        double diagonal_value = A;
        struct edge *ptrl = ptr->list_edges;

        while (ptrl != NULL)
        {
            uint32_t v = ptrl->id;
            uint32_t link_type = ptrl->link_type;

            // Neighbour volume is link by the citoplasm
            if (link_type == 0)
            {
                value = -B;
                diagonal_value += B;
            }
            // Neighbour volume is link by a gap junction
            else if (link_type == 1)
            {
                value = -C;
                diagonal_value += C;
            }
            coeff.push_back(Eigen::Triplet<double>(u,v,value));

            ptrl = ptrl->next;
        }
        coeff.push_back(Eigen::Triplet<double>(u,u,diagonal_value));
        ptr = ptr->next;
    }

    //for (int i = 0; i < coeff.size(); i++)
    //    printf("(%d,%d) = %.10lf\n",coeff[i].row(),coeff[i].col(),coeff[i].value());
    a.setFromTriplets(coeff.begin(),coeff.end());
    a.makeCompressed();

    return a;
}

void solve_monodomain (struct monodomain_solver *the_monodomain_solver,\
                        struct ode_solver *the_ode_solver,\
                        struct graph *the_purkinje_network,\
                        struct user_options *options)
{
    assert(the_monodomain_solver);
    assert(the_ode_solver);
    assert(the_purkinje_network);
    assert(options);

    set_ode_initial_condition_for_all_volumes (the_ode_solver);
    //print_state_vector(the_ode_solver);

    // Assemble the matrix
    uint32_t nc = the_ode_solver->n_active_cells;

    Eigen::SparseMatrix<double> A(nc,nc);
    A = assemble_matrix(the_monodomain_solver,the_purkinje_network);
    Eigen::SparseLU< Eigen::SparseMatrix<double> > sparse_solver(A);

    // Declare RHS and the solution vector
    Eigen::VectorXd b(nc);
    Eigen::VectorXd x(nc);

    double tmax = the_monodomain_solver->tmax;
    double dt = the_monodomain_solver->dt;
    uint32_t M = tmax / dt;

    // Time loop
    for (int i = 0; i < M; i++)
    {
        double t = i*dt;

        // Write the solution to .vtk file
        if (i % 100 == 0) 
            write_to_VTK(the_purkinje_network,the_ode_solver,i);
        
        // Solve the PDE (diffusion phase)
        //assembleLoadVector(b);
        //x = sparseSolver.solve(b);
        //moveVstar(x);

        // Solve the ODEs (reaction phase)
        //solveODE(t);

        // Jump to the next iteration
        //nextTimestep();
    }

}

void write_to_VTK (struct graph *the_purkinje_network, struct ode_solver *the_ode_solver, int iter)
{
    char filename[50];
    FILE *file;

    uint32_t n_odes = the_ode_solver->num_ode_equations;
    uint32_t np = the_purkinje_network->total_nodes;
    uint32_t ne = the_purkinje_network->total_edges;
    struct node *ptr = the_purkinje_network->list_nodes;
    double *sv = the_ode_solver->sv;

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
    ptr = the_purkinje_network->list_nodes;
    while (ptr != NULL)
    {
        struct edge *ptrl = ptr->list_edges;
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

    for (int i = 0; i < np; i++)
        fprintf(file,"%.10lf\n",sv[i*n_odes]);
    
    fclose(file);
}