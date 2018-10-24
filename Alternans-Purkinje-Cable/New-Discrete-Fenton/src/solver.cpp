#include "../include/solver.h"

Solver::Solver (Options *user_options)
{
    nthreads = user_options->num_threads;
    dt = user_options->dt;
    tmax = user_options->tmax;
    print_rate = user_options->print_rate;
    sst_rate = user_options->sst_rate;
    network_filename = user_options->network_filename;
    steady_filename = user_options->sst_filename;
    plot_filename = user_options->plot_filename;
    diameter = user_options->start_diameter;
    sigma_c = user_options->sigma_c;
    g_GAP = user_options->G_gap;

    dx = user_options->start_h / user_options->num_div_cell;
    M = nearbyint(tmax/dt);
    beta = 4.0 / diameter * 1.0e-04;

    stim_config = new Stimulus(user_options);
    //stim_config->print();

    the_purkinje_network = set_purkinje_mesh_from_file(network_filename,dx);
    //the_purkinje_network->print();

    the_purkinje_network->set_gap_junctions(user_options->num_div_cell);
    //the_purkinje_network->print();

    set_plot_cells();
    set_celular_model(stim_config);
    set_control_volumes();
    set_derivative();
    set_functions();
    set_velocity_points();
    set_plot_points();

    if (user_options->steady_state)
        set_initial_conditions_from_file();
    else
        set_initial_conditions_from_default();
    
}

void Solver::solve ()
{
    // Open the SST that will store the steady_state variables of the system on the time SST_RATE
    FILE *sst_file = fopen("output.sst","w+");

    int np = the_purkinje_network->get_total_nodes();
    
    // Build the discretization matrix
    SpMat A(np,np);
    set_matrix(A);
    // Apply a LU decomposition over the matrix
    SparseLU<SpMat> sparseSolver(A);

    // Declare RHS and the solution vector
    VectorXd b(np);
    VectorXd x(np);
    
    // Time loop
    for (int i = 0; i < M; i++)
    {
        double t = i*dt;

        for (int k = 0; k < 100; k++)

        // Print the progress of the solution
        #ifdef OUTPUT
        print_progress(i,M);
        #endif

        // Write the solution data to a file
        //writePlotData(t);
        
        //if (i % 10 == 0)
        //    writeStateVectorToFile(i,max_num_digits,t);

        // Write the solution to .vtk file
        #ifdef VTK
        //if (i % PRINT_RATE == 0) writeVTKFile(i);
        //if (i == SST_RATE-1) writeSteadyStateFile(sstFile);
        #endif

        // Solve the PDE (diffusion phase)
        //assembleLoadVector(b);
        //x = sparseSolver.solve(b);

        //moveVstar(x);

        // Solve the ODEs (reaction phase)
        //solveODE(t);

        // Calculate the maximum derivative for each volume
        //calcMaxDerivative(t);

        // Jump to the next iteration
        //nextTimestep();
    }


    fclose(sst_file);
}

// Build the coefficient matrix considering cells with the same diameter
void Solver::set_matrix (SpMat &a)
{
    cout << "[Solver] Building discretization matrix" << endl;
    // Compute the coefficients values
    double A = (4.0*g_GAP*dx) / (M_PI*diameter*diameter);
    double B = (sigma_c);
    double C = (beta*Cm*dx*dx) / (dt);

    //printf("A = %.10lf\n",A);
    //printf("B = %.10lf\n",B);
    //printf("C = %.10lf\n",C);

    // Non-zero coefficients
    vector<T> coeff;

    double diagonal_value;
    Node *ptr = the_purkinje_network->get_list_nodes();
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
            coeff.push_back(T(u,v,value));

            ptrl = ptrl->next;
        }
        coeff.push_back(T(u,u,diagonal_value));

        ptr = ptr->next;
    }
/*
    FILE *file = fopen("matrix.txt","w+");
    fprintf(file,"%d %d\n",g->get_total_nodes(),g->get_total_nodes());
    for (size_t i = 0; i < coeff.size(); i++)
        fprintf(file,"%d %d %.10lf\n",coeff[i].row(),coeff[i].col(),coeff[i].value());
    fclose(file);
*/
    a.setFromTriplets(coeff.begin(),coeff.end());
    a.makeCompressed();

}

void Solver::set_plot_cells ()
{
    // Allocate and initialize the vector 'ids' with the volumes to be plotted from the input file .plt
    plot = (Plot*)malloc(sizeof(Plot));
    FILE *pltFile = fopen(plot_filename.c_str(),"r");
    if (!pltFile)
    {
        cerr << "[-] ERROR! Cannot open PLT file!" << endl;
        exit(EXIT_FAILURE);
    }

    if (!fscanf(pltFile,"%d",&plot->np)) 
    {
        cerr << "[-] ERROR! Reading PLT file!" << endl;
        exit(EXIT_FAILURE);
    }
    plot->ids = (int*)malloc(sizeof(int)*plot->np);

    for (int i = 0; i < plot->np; i++)
    {
        if (!fscanf(pltFile,"%d",&plot->ids[i]))
        {
            cerr << "[-] ERROR! Reading PLT file!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    fclose(pltFile);
}

void Solver::set_control_volumes ()
{
    // Capture the number of equations of the celullar model
    Node *ptr = the_purkinje_network->get_list_nodes();
    int neq = NUM_EQ;
    int np = the_purkinje_network->get_total_nodes();
    vol = (CVolume*)malloc(sizeof(CVolume)*np);
    for (int i = 0; i < np; i++, ptr = ptr->next)
    {
        vol[i].type = ptr->type;
        vol[i].yOld = (double*)calloc(neq,sizeof(double));
        vol[i].yNew = (double*)calloc(neq,sizeof(double));
        vol[i].yStar = (double*)calloc(neq,sizeof(double));
    }
}

void Solver::set_derivative ()
{
    int np = the_purkinje_network->get_total_nodes();
    dvdt = (Derivative*)malloc(np*sizeof(Derivative));
    for (int i = 0; i < np; i++) 
        dvdt[i].value = 0;
}

void Solver::set_functions ()
{
    func = (Func*)malloc(sizeof(Func)*NUM_EQ);  
    func[0] = dvdt__Nob;
    func[1] = dmdt__Nob;
    func[2] = dhdt__Nob;
    func[3] = dndt__Nob;
}

void Solver::set_initial_conditions_from_file ()
{
    cout << "[Solver] Setting initial conditions from SST file: " << steady_filename << endl;
    
    FILE *sstFile = fopen(steady_filename.c_str(),"r");
    if (!sstFile)
    {
        cerr << "[-] ERROR! Could not open SST file!" << endl;
        exit(EXIT_FAILURE);
    }

    int neq = NUM_EQ;
    int np = the_purkinje_network->get_total_nodes();
    // Iterate over all the points
    for (int i = 0; i < np; i++)
    {
        // Iterate over all the ODEs equations
        for (int j = 0; j < neq; j++)
            if (!fscanf(sstFile,"%lf",&vol[i].yOld[j]))
            {
                cerr << "[-] ERROR! Could not open SST file!" << endl;
                exit(EXIT_FAILURE); 
            } 
    }
    fclose(sstFile);
}

void Solver::set_initial_conditions_from_default ()
{
    cout << "[Solver] Setting default initial conditions" << endl;

    int neq = NUM_EQ;
    int np = the_purkinje_network->get_total_nodes();
    for (int i = 0; i < np; i++)
        for (int j = 0; j < neq; j++)
            vol[i].yOld[j] = YO__NOBLE[j];
}

void Solver::set_velocity_points ()
{
    vel = (Velocity*)malloc(sizeof(Velocity));
    vel->velocityFile = fopen("output/velocity.txt","w+");

    // First point is always the source
    vel->np = plot->np-1;
    vel->id_source = plot->ids[0];
    vel->ids = (int*)malloc(sizeof(int)*vel->np);
    vel->t2 = (double*)malloc(sizeof(double)*vel->np);
    for (int i = 0; i < vel->np; i++) 
        vel->ids[i] = plot->ids[i+1];
}

void Solver::set_plot_points ()
{
    char filename[200];
    plot->plotFile = (FILE**)malloc(sizeof(FILE*)*(plot->np-1));
    for (int i = 1; i < plot->np; i++)
    {
        plot->plotFile[i-1] = (FILE*)malloc(sizeof(FILE));
        sprintf(filename,"output/data-%d.dat",plot->ids[i]);
        plot->plotFile[i-1] = fopen(filename,"w+");
    }
}

void Solver::print ()
{
    cout << "////////////////////////////////////////////////////////////////" << endl;
    cout << "number of threads = " << nthreads << endl;
    cout << "dt = " << dt << endl;
    cout << "tmax = " << tmax << endl;
    cout << "number of timesteps = " << M << endl;
    cout << "dx = " << dx << endl;
    cout << "network_filename = " << network_filename << endl;
    cout << "steady_state_filename = " << steady_filename << endl;
    cout << "plot_filename = " << plot_filename << endl;
    cout << "print_rate = " << print_rate << endl;
    cout << "sst_rate = " << sst_rate << endl;
    cout << "diameter = " << diameter << endl;
    cout << "beta = " << beta << endl;
    cout << "Cm = " << Cm << endl;
    cout << "sigma_c = " << sigma_c << endl;
    cout << "G_gap = " << g_GAP << endl;
    cout << "////////////////////////////////////////////////////////////////" << endl;
}

void print_progress (int iter, int max_iter)
{
    double progress = iter / (double)max_iter;
    
    cout << "Progress: " << int(progress * 100.0) << " %\r";
    cout.flush();
}