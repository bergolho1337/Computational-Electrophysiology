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
    
    cout << "[Solver] Time loop" << endl;
    
    // Time loop
    for (int i = 0; i < M; i++)
    {
        double t = i*dt;

        // Print the progress of the solution
        #ifdef OUTPUT
        print_progress(i,M);
        #endif

        // Write the solution data to a file
        write_plot_data(t);
        
        //if (i % 10 == 0)
        //    writeStateVectorToFile(i,max_num_digits,t);

        // Write the solution to .vtk file
        #ifdef VTK
        if (i % print_rate == 0) write_VTK_file(i);
        if (i == sst_rate-1) write_steady_state_file(sst_file);
        #endif

        // Solve the PDE (diffusion phase)
        assemble_load_vector(b);
        x = sparseSolver.solve(b);

        move_Vstar(x);

        // Solve the ODEs (reaction phase)
        solve_ODE(t);

        // Calculate the maximum derivative for each volume
        calc_max_derivative(t,stim_config->start_period);

        // Jump to the next iteration
        next_timestep();
    }

    calc_velocity();

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

void Solver::assemble_load_vector (VectorXd &b)
{

    double C = (beta*Cm*dx*dx) / (dt);

    int np = b.size();
    for (int i = 0; i < np; i++)
        b(i) = vol[i].yOld[0] * C;
    
    /*
    FILE *file = fopen("rhs.txt","w+");
    fprintf(file,"%d\n",b.size());
    for (size_t i = 0; i < b.size(); i++)
        fprintf(file,"%.10lf\n",b(i));
    fclose(file);
    */

}

void Solver::move_Vstar (const VectorXd vm)
{
    int neq = NUM_EQ;
    int np = vm.size();
    for (int i = 0; i < np; i++)
    {
        vol[i].yStar[0] = vm(i);
        for (int j = 1; j < neq; j++)
            vol[i].yStar[j] = vol[i].yOld[j];
    }
}

void Solver::solve_ODE (double t)
{
    int neq = NUM_EQ;
    int ncells = the_purkinje_network->get_total_nodes();
    #pragma omp parallel for num_threads(nthreads)
    for (int id = 0; id < ncells; id++)
    {
        // V^n+1 = V^n+1/2 + f*dt
        double f = func[0](vol[id].type,id,t,vol[id].yStar);
        vol[id].yNew[0] = vol[id].yStar[0] + f*dt;
        // gate^n+1 = gate^n + dt*f
        for (int j = 1; j < neq; j++)
        {
            f = func[j](vol[id].type,id,t,vol[id].yOld);
            vol[id].yNew[j] = vol[id].yOld[j] + f*dt;
        } 
    }
}

void Solver::calc_max_derivative (double t, double current_period)
{
    int np = the_purkinje_network->get_total_nodes();
    for (int i = 0; i < np; i++)
    {
        double diff = vol[i].yNew[0] - vol[i].yOld[0];
        // Considering the calculus after the first stimulus
        if (diff > dvdt[i].value && t > 0.0 && t < current_period)
        {
            dvdt[i].value = diff;
            dvdt[i].t = t;
        }
    }
}

void Solver::calc_velocity ()
{
    cout << "[Solver] Calculating propagation velocity" << endl;

    FILE *vFile = fopen("output/v.txt","w+");
    int np = vel->np;
    double *dist = the_purkinje_network->get_dist();
    
    for (int i = 0; i < np; i++)
    {    
        // Compute the distance from the source
        the_purkinje_network->dijkstra(vel->ids[i]);
        dist = the_purkinje_network->get_dist();
        double t = dvdt[vel->ids[i]].t - dvdt[vel->ids[i] - OFFSET].t;
        double velocity = dist[vel->ids[i] - OFFSET] / t;

        // Check if the value is beyond the tolerance
	    if (t < 0 || fabs(dvdt[vel->ids[i]].value - dvdt[vel->ids[i] - OFFSET].value) > 16.0)
            velocity = 0.0;
        fprintf(vel->velocityFile,"\n\n[!] Propagation velocity! Id = %d\n",vel->ids[i]);
        fprintf(vel->velocityFile,"t1 = %.10lf\n",dvdt[vel->ids[i] - OFFSET].t);
        fprintf(vel->velocityFile,"dvdt[%d] = %.10lf\n\n",vel->ids[i]- OFFSET ,dvdt[vel->ids[i] - OFFSET].value);
        fprintf(vel->velocityFile,"t2 = %.10lf\n",dvdt[vel->ids[i]].t);
        fprintf(vel->velocityFile,"dvdt[%d] = %.10lf\n",vel->ids[i],dvdt[vel->ids[i]].value);
        fprintf(vel->velocityFile,"delta_x = %.10lf\n",dist[vel->ids[i] - OFFSET]);
        fprintf(vel->velocityFile,"delta_t = %.10lf\n",t);
        fprintf(vel->velocityFile,"\n!!!!!!!! Propagation velocity = %lf cm/s !!!!!!!!!!\n",velocity*1000.0);
        fprintf(vel->velocityFile,"\n=============================================================\n\n");

        fprintf(vFile,"%lf\n",velocity*1000.0);
    }
    fclose(vel->velocityFile);
    fclose(vFile);
}

void Solver::next_timestep ()
{
    int np = the_purkinje_network->get_total_nodes();
    for (int i = 0; i < np; i++) 
        swap(&vol[i].yOld,&vol[i].yNew);
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

void Solver::write_plot_data (double t)
{
    for (int i = 1; i < plot->np; i++)
        fprintf(plot->plotFile[i-1],"%.10lf %.10lf\n",t,vol[plot->ids[i]].yOld[0]);
}

void Solver::write_VTK_file (int iter)
{
    FILE *file;
    int np, ne;
    char filename[50];
    Node *ptr = the_purkinje_network->get_list_nodes();
    np = the_purkinje_network->get_total_nodes();
    ne = the_purkinje_network->get_total_edges();

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
    ptr = the_purkinje_network->get_list_nodes();
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
    ptr = the_purkinje_network->get_list_nodes();
    while (ptr != NULL)
    {
        fprintf(file,"%e\n",vol[ptr->id].yOld[0]);
        ptr = ptr->next;
    }
    fclose(file);
}

void Solver::write_steady_state_file (FILE *sstFile)
{
    int neq = NUM_EQ;
    int np = the_purkinje_network->get_total_nodes();
    for (int i = 0; i < np; i++)
    {
        for (int j = 0; j < neq; j++)
            fprintf(sstFile,"%.10lf ",vol[i].yOld[j]);
        fprintf(sstFile,"\n");
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

void swap (double **a, double **b)
{
    double *tmp = *a;
    *a = *b;
    *b = tmp;
}