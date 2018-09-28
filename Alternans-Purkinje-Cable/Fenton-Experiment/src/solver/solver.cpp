#include "solver.h"

static const char *opt_string =   "c:o:n:p:t:e:f:h";

static const struct option long_options[] = {
        { "config_file", required_argument, NULL, 'c' },
        { "output_dir", required_argument , NULL, 'o' },
        { "num_threads", required_argument, NULL, 'n' },
        { "print_rate",required_argument , NULL, 'p' },
        { "dt", required_argument, NULL, 'e' },
        { "simulation_time", required_argument, NULL, 'f' },
        { "sigma", required_argument, NULL, SIGMA},
        { "beta", required_argument, NULL, BETA},
        { "cm", required_argument, NULL, CM},
        { "stimulus", required_argument, NULL, STIM_OPT},
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
};

// Number of threads to solve the system of ODEs
//static constexpr int nthreads = 2;
/*
Solver::Solver (int argc, char *argv[])
{
    dt = atof(argv[2]);
    tmax = atof(argv[3]);
    mesh_filename = argv[4];
    steady_filename = argv[5];
    plot_filename = argv[6];
    setTypeCell();
    M = nearbyint(tmax/dt);
    g = new Graph(mesh_filename,dx);
    //g->dijkstra(plot->ids[0]);

    setControlVolumes();
    setDerivative();
    setFunctions();
    setInitCondFromFile();
    setVelocityPoints();
    setPlot(); 
    setSensibilityParam(argc,argv);

    //g->print();   
}

void Solver::solve ()
{
    #ifdef OUTPUT
    printf("[!] Solving transient problem ... \n");
    printf("[!] Progress\n");
    fflush(stdout);
    #endif

    int np = g->getTotalNodes();

    // Build the matrix
    SpMat A(np,np);
    #ifdef DIAMETER
    setMatrix2(A);
    #else
    setMatrix(A);
    #endif
    SparseLU<SpMat> sparseSolver(A);
    
    // Declare RHS and the solution vector
    VectorXd b(np);
    VectorXd x(np);
    
    // Time loop
    for (int i = 0; i < M; i++)
    {
        double t = i*dt;

        // Print the progress of the solution
        #ifdef OUTPUT
        printProgress(i,M);
        #endif

        // Write the solution data to a file
        writePlotData(t);

        // Write the solution to .vtk file
        #ifdef VTK
        if (i % 10 == 0) writeVTKFile(i);
        #endif

        // Solve the PDE (diffusion phase)
        assembleLoadVector(b);
        x = sparseSolver.solve(b);
        moveVstar(x);

        // Solve the ODEs (reaction phase)
        solveODE(t);

        // Calculate the maximum derivative for each volume
        calcMaxDerivative(t);

        // Jump to the next iteration
        nextTimestep();
    }
    #ifdef OUTPUT
    printf("ok\n");
    #endif

    // Calculate the propagation velocity for the plot ids
    calcVelocity();
}

void Solver::setSensibilityParam (int argc, char *argv[])
{
    // Default values
    if (argc-1 == 6)
    {
        alfa = 1.375;
        d1 = 0.002;
        SIGMA = 0.004;
    }
    else if (argc-1 == 8)
    {
        alfa = atof(argv[7]);
        d1 = atof(argv[8]);
        SIGMA = 0.004;
    }
    // User-defined
    else
    {
        alfa = atof(argv[7]);
        d1 = atof(argv[8]);
        SIGMA = atof(argv[9]);
    }
    BETA = 4.0 / d1 * 1.0e-04;
}

void Solver::setControlVolumes ()
{
    // Capture the number of equations of the celullar model
    Node *ptr = g->getListNodes();
    int neq = num_eq;
    int np = g->getTotalNodes();
    vol = (CVolume*)malloc(sizeof(CVolume)*np);
    for (int i = 0; i < np; i++, ptr = ptr->next)
    {
        vol[i].type = ptr->type;
        vol[i].yOld = (double*)calloc(neq,sizeof(double));
        vol[i].yNew = (double*)calloc(neq,sizeof(double));
        vol[i].yStar = (double*)calloc(neq,sizeof(double));
    }
}

void Solver::setFunctions ()
{
  func = (Func*)malloc(sizeof(Func)*num_eq);  
  func[0] = dvdt__Nob;
  func[1] = dmdt__Nob;
  func[2] = dhdt__Nob;
  func[3] = dndt__Nob;
}

void Solver::setTypeCell ()
{
    // Allocate and initialize the vector 'ids' with the volumes to be plotted from the input file .plt
    plot = (Plot*)malloc(sizeof(Plot));
    FILE *pltFile = fopen(plot_filename.c_str(),"r");

    if (!fscanf(pltFile,"%d",&plot->np)) error("Reading PLT file");
    plot->ids = (int*)malloc(sizeof(int)*plot->np);

    for (int i = 0; i < plot->np; i++)
        if (!fscanf(pltFile,"%d",&plot->ids[i])) error("Reading PLT file");

    fclose(pltFile);
}

void Solver::setInitCondFromFile ()
{
    FILE *sstFile = fopen(steady_filename.c_str(),"r");
    if (!sstFile) error("Could not open SST file");
    int neq = num_eq;
    int np = g->getTotalNodes();
    // Iterate over all the points
    for (int i = 0; i < np; i++)
    {
        // Iterate over all the ODEs equations
        for (int j = 0; j < neq; j++)
            if (!fscanf(sstFile,"%lf",&vol[i].yOld[j])) error("Reading SST file.");
    }
    fclose(sstFile);
}

void Solver::setVelocityPoints ()
{
    vel = (Velocity*)malloc(sizeof(Velocity));
    vel->velocityFile = fopen("Output/velocity.txt","w+");
    // First point is always the source
    vel->np = plot->np-1;
    vel->id_source = plot->ids[0];
    vel->ids = (int*)malloc(sizeof(int)*vel->np);
    vel->t2 = (double*)malloc(sizeof(double)*vel->np);
    for (int i = 0; i < vel->np; i++) 
        vel->ids[i] = plot->ids[i+1];
}

void Solver::setPlot ()
{
    char filename[200];
    plot->plotFile = (FILE**)malloc(sizeof(FILE*)*(plot->np-1));
    for (int i = 1; i < plot->np; i++)
    {
        plot->plotFile[i-1] = (FILE*)malloc(sizeof(FILE));
        sprintf(filename,"Output/data%d.dat",plot->ids[i]);
        plot->plotFile[i-1] = fopen(filename,"w+");
    }
}

void Solver::setDerivative ()
{
    int np = g->getTotalNodes();
    dvdt = (Derivative*)malloc(np*sizeof(Derivative));
    for (int i = 0; i < np; i++) dvdt[i].value = 0;
}

void Solver::setMatrix (SpMat &a)
{
    // Compute the coefficients values
    double A = 4.0 / (RPMJ*M_PI*d1*d1*dx);
    double B = (SIGMA) / (dx*dx);
    double C = (BETA*Cm) / (dt);
    double D = (BETA*Cm*alfa) / (dt);
    double E = (BETA*Cm*dx*dx) / (dt);

    // Non-zero coefficients
    vector<T> coeff;

    Node *ptr = g->getListNodes();
    while (ptr != NULL)
    {
        int u = ptr->id;
        int type = ptr->type;
        Edge *ptrl = ptr->edges;
        
        // PMJ
        if (type == 1)
        {
            double value = -1.0 / D;
            while (ptrl != NULL)
            {
                int v = ptrl->dest->id;
                coeff.push_back(T(u,v,value));
                ptrl = ptrl->next;
            }
            value = (1.0 + D) / D;
            coeff.push_back(T(u,u,value)); 
        }
        // Purkinje cell
        else
        {
            bool isPMJ = isConnToPMJ(ptr->edges);
            //Not link to a PMJ, so normal edge with a Purkinje cell
            if (isPMJ == false)
            {
                double value = -SIGMA / E;
                while (ptrl != NULL)
                {
                    int v = ptrl->dest->id;
                    coeff.push_back(T(u,v,value));
                    ptrl = ptrl->next;
                }
                value = (ptr->num_edges*SIGMA + E) / E;
                coeff.push_back(T(u,u,value));
            }
            // Is a special link to a Purkinje cell - PMJ
            else
            {
                double sum = C;
                while (ptrl != NULL)
                {
                    int v = ptrl->dest->id;
                    if (ptrl->dest->type == 0)
                    {
                        double value = -B / C;
                        sum += B;
                        coeff.push_back(T(u,v,value));
                    }
                    else
                    {
                        double value = -A / C;
                        sum += A;
                        coeff.push_back(T(u,v,value));
                    }
                    ptrl = ptrl->next;
                }
                sum /= C;
                coeff.push_back(T(u,u,sum));
            }  
        }
        ptr = ptr->next;
    }
    a.setFromTriplets(coeff.begin(),coeff.end());
    a.makeCompressed();
}

void Solver::setMatrix2 (SpMat &a)
{
    // Compute the coefficients values
    double B = (SIGMA) / (dx*dx);
    
    // Non-zero coefficients
    vector<T> coeff;

    Node *ptr = g->getListNodes();
    while (ptr != NULL)
    {
        int u = ptr->id;
        int type = ptr->type;
        double d_u = ptr->d;
        BETA = 4.0 / d_u * 1.0e-04;
        Edge *ptrl = ptr->edges;
        
        // PMJ
        if (type == 1)
        {
            double D = (BETA*Cm*alfa) / (dt);
            double value = -1.0 / D;
            while (ptrl != NULL)
            {
                int v = ptrl->dest->id;
                coeff.push_back(T(u,v,value));
                ptrl = ptrl->next;
            }
            value = (1.0 + D) / D;
            coeff.push_back(T(u,u,value)); 
        }
        // Purkinje cell
        else
        {
            bool isPMJ = isConnToPMJ(ptr->edges);
            //Not link to a PMJ, so normal edge with a Purkinje cell
            if (isPMJ == false)
            {
                double E = (BETA*Cm*d_u*d_u*dx*dx) / (dt);
                double sum = 1.0;
                while (ptrl != NULL)
                {
                    int v = ptrl->dest->id;
                    double d_v = ptrl->dest->d;
                    double value = (SIGMA*d_v*d_v) / E;
                    coeff.push_back(T(u,v,-value));
                    sum += value;
                    ptrl = ptrl->next;
                }
                coeff.push_back(T(u,u,sum));
            }
            // Is a special link to a Purkinje cell - PMJ
            else
            {
                double C = (BETA*Cm) / (dt);
                double sum = C;
                while (ptrl != NULL)
                {
                    int v = ptrl->dest->id;
                    double d_v = ptrl->dest->d;
                    // Purkinje cell - Purkinje cell
                    if (ptrl->dest->type == 0)
                    {
                        double value = -B / C;
                        sum += B;
                        coeff.push_back(T(u,v,value));
                    }
                    // Purkinje cell - PMJ
                    else
                    {
                        double A = (4.0) / (RPMJ*M_PI*d_v*d_v*dx);
                        double value = A / C;
                        sum += A;
                        coeff.push_back(T(u,v,-value));
                    }
                    ptrl = ptrl->next;
                }
                sum /= C;
                coeff.push_back(T(u,u,sum));
            }  
        }
        ptr = ptr->next;
    }
    a.setFromTriplets(coeff.begin(),coeff.end());
    a.makeCompressed();
}

void Solver::assembleLoadVector (VectorXd &b)
{
    int np = b.size();
    for (int i = 0; i < np; i++)
        b(i) = vol[i].yOld[0];
}

void Solver::moveVstar (const VectorXd vm)
{
    int neq = num_eq;
    int np = vm.size();
    for (int i = 0; i < np; i++)
    {
        vol[i].yStar[0] = vm(i);
        for (int j = 1; j < neq; j++)
            vol[i].yStar[j] = vol[i].yOld[j];
    }
}

void Solver::solveODE (double t)
{
    int neq = num_eq;
    int ncells = g->getTotalNodes();
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

void Solver::writePlotData (double t)
{
    for (int i = 1; i < plot->np; i++)
        fprintf(plot->plotFile[i-1],"%.10lf %.10lf\n",t,vol[plot->ids[i]].yOld[0]);
}

void Solver::writeVTKFile (int iter)
{
    FILE *file;
    int np, ne;
    char filename[50];
    Node *ptr = g->getListNodes();
    np = g->getTotalNodes();
    ne = g->getTotalEdges();

    // Write the transmembrane potential
    sprintf(filename,"VTK/sol%d.vtk",iter);
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
    ptr = g->getListNodes();
    while (ptr != NULL)
    {
        Edge *ptrl = ptr->edges;
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
    ptr = g->getListNodes();
    while (ptr != NULL)
    {
        fprintf(file,"%e\n",vol[ptr->id].yOld[0]);
        ptr = ptr->next;
    }
    fclose(file);
}

bool Solver::isConnToPMJ (Edge *e)
{
    while (e != NULL)
    {
        if (e->dest->type == 1) return true;
        e = e->next;
    }
    return false;
}

void Solver::calcMaxDerivative (double t)
{
    int np = g->getTotalNodes();
    for (int i = 0; i < np; i++)
    {
        double diff = vol[i].yNew[0] - vol[i].yOld[0];
        // Considering the calculus after the first stimulus
        if (diff > dvdt[i].value && t > 0.0 && t < 250.0)
        {
            dvdt[i].value = diff;
            dvdt[i].t = t;
        }
    }
}

void Solver::calcVelocity ()
{
    FILE *vFile = fopen("Output/v.txt","w+");
    FILE *dFile = fopen("Output/delay.txt","w+");
    FILE *tFile = fopen("Output/time.txt","w+");
    int np = vel->np;
    int nterm = g->getNTerm();
    int *term = g->getTerm();

    g->dijkstra(0);
    double *dist = g->getDist();
    for (int i = 0; i < np; i++)
    {
        printf("Dist = %.10lf\n",dist[vel->ids[i]]);
    }
    for (int i = 0; i < np; i++)
    {    
        // Compute the distance from the source
        g->dijkstra(vel->ids[i]);
        dist = g->getDist();
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
        fprintf(tFile,"%lf\n",dvdt[vel->ids[i]].t);
    }
    for (int i = 0; i < nterm; i++)
    { 
        // Calculate the delay between the terminal and 10 volumes behind it
        printf("Delay between %d and %d\n",term[i],term[i]-10);
        fprintf(dFile,"%.10lf\n",dvdt[term[i]].t-dvdt[term[i]-10].t);
    }
    //fprintf(dFile,"%.10lf\n",dvdt[plot->delay[1]].t-dvdt[plot->ids[0]].t);
    fclose(vel->velocityFile);
    fclose(dFile);
    fclose(vFile);
    fclose(tFile);
}

void Solver::nextTimestep ()
{
    int np = g->getTotalNodes();
    for (int i = 0; i < np; i++) swap(&vol[i].yOld,&vol[i].yNew);
}

void Solver::error (const char msg[])
{
    printf("[-] ERROR ! %s !\n",msg);
    exit (EXIT_FAILURE);
}
*/

void print_user_options (struct user_options *options)
{
    cout << PRINT_LINE << endl;
    cout << "num_threads = " << options->num_threads << endl;
    cout << "final_time = " << options->final_time << endl;
    cout << "dt = " << options->dt << endl;
    cout << "output_directory = " << options->out_dir_name << endl;
    cout << "print_rate = " << options->print_rate << endl;
    cout << "output_steady_state_directory = " << options->out_steady_state_dir << endl;
    cout << "steady_state_filename = " << options->steady_state_filename << endl;
    cout << "sigma = " << options->sigma << endl;
    cout << "cm = " << options->cm << endl;
    cout << "beta = " << options->beta << endl;
    cout << "config_file = " << options->config_file << endl;
    cout << PRINT_LINE << endl;
    cout << "name = " << options->purkinje_config->name << endl;
    cout << "network_file = " << options->purkinje_config->network_filename << endl;
    cout << "start_discretization = " << options->purkinje_config->start_h << endl;
    cout << PRINT_LINE << endl;
    cout << "stim_start = " << options->stim_configs->stim_start << endl;
    cout << "stim_current = " << options->stim_configs->stim_current << endl;
    cout << "stim_duration = " << options->stim_configs->stim_duration << endl;
    cout << "period_step = " << options->stim_configs->period_step << endl;
    cout << "start_period = " << options->stim_configs->start_period << endl;
    cout << "end_period = " << options->stim_configs->end_period << endl;
    cout << "n_cycles = " << options->stim_configs->n_cycles << endl;
    cout << "id_limit = " << options->stim_configs->id_limit << endl;
    cout << PRINT_LINE << endl;
}

void usage (const char pName[])
{
    cerr << PRINT_LINE_2 << endl;
    cerr << "Usage:> " << pName << " -c <configuration_file>" << endl; 
    cerr << PRINT_LINE << endl;
    cerr << "<configuration_file> = Configuration file with the parameters of the simulation (.ini)" << endl;
    cerr << PRINT_LINE << endl;
    cerr << PRINT_LINE_2 << endl;
    /*
     // Old Usage
    fprintf(stderr,"=====================================================================================\n");
    fprintf(stderr,"Usage:> %s -t/s <dt> <tmax> <mesh_file> <steady_state_file> <plot_file> [ALPHA] [d1] [SIGMA]\n",pName);
    fprintf(stderr,"-------------------------------------------------------------------------------------\n");
    fprintf(stderr,"-t = Steady-State mode\n");
    fprintf(stderr,"-s = Solver mode\n\n");
    fprintf(stderr,"<dt> = Size of the time discretization\n");
    fprintf(stderr,"<tmax> = Maximum simulation time\n");
    fprintf(stderr,"<mesh_file> = File with the mesh points and conections\n");
    fprintf(stderr,"<steady_state_file> = File with the steady state solution\n");
    fprintf(stderr,"<plot_file> = File with ids to be plotted\n");
    fprintf(stderr,"-------------------------------------------------------------------------------------\n");
    fprintf(stderr," !! Optional parameters !! (For the sensibility analysys)\n");
    fprintf(stderr,"[ALPHA] = R_pmj*Vol_pmj\n");
    fprintf(stderr,"[d1] = Diameter Purkinje cell\n");
    fprintf(stderr,"[SIGMA] = Conductivity Gap Junction + Citoplasm\n");
    fprintf(stderr,"-------------------------------------------------------------------------------------\n");
    fprintf(stderr,"Steady-State Example: %s -t 0.1 5000 cable_dog.msh steady_dog.sst cable_dog.plt 1.375 0.002 0.004\n",pName);
    fprintf(stderr,"Solver Example: %s -s 0.1 1000 cable_dog.msh steady_dog.sst cable_dog.plt 1.375 0.002 0.004\n",pName);
    fprintf(stderr,"-------------------------------------------------------------------------------------\n");
    */
}

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *options = (struct user_options*)malloc(sizeof(struct user_options));

    options->num_threads = 1;
    options->num_threads_was_set = false;

    options->out_dir_name = NULL;
    options->out_dir_name_was_set = false;

    options->final_time = 10.0;
    options->final_time_was_set = false;

    options->print_rate = 1;
    options->print_rate_was_set = false;

    options->out_steady_state_dir = NULL;
    options->out_steady_state_dir_was_set = false;

    options->dt = 0.01;
    options->dt_was_set = false;

    options->steady_state_filename = NULL;
    options->steady_state_filename_was_set = false;

    options->sigma = 0.0000176;
    options->sigma_was_set = false;

    options->cm = 1.2;
    options->cm_was_set = false;

    options->beta = 0.14;
    options->beta_was_set = false;

    options->config_file = NULL;

    options->purkinje_config = NULL;
    options->stim_configs = NULL;

    return options;
}

void solve_model (int argc, char *argv[])
{
    struct user_options *options = new_user_options(argc,argv);
    get_config_file (argc,argv,options);

    if (options->config_file) 
    {
        // Here we parse the config file
        if (ini_parse (options->config_file, parse_config_file, options) < 0) 
        {
            fprintf (stderr, "Error: Can't load the config file %s\n", options->config_file);
            exit(EXIT_FAILURE);
        }
    }

    print_user_options(options);

}

void get_config_file (int argc, char *argv[], struct user_options *user_args)
{
    int opt = 0;

    int option_index;

    opt = getopt_long (argc, argv, opt_string, long_options, &option_index);

    while (opt != -1) 
    {
        switch (opt) 
        {
            case 'c':
                user_args->config_file = optarg;
                return;
            default:
                break;
        }
        opt = getopt_long (argc, argv, opt_string, long_options, &option_index);
    }

    // We reset the index after parsing the config_file
    optind = 1;
}

int parse_config_file (void *user, const char *section, const char *name, const char *value) 
{
    struct user_options *pconfig = (struct user_options*) user;

    if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "num_threads")) 
    {
        pconfig->num_threads = (int)strtol (value, NULL, 10);
        pconfig->num_threads_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "dt")) 
    {
        pconfig->dt = strtod(value, NULL);
        pconfig->dt_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "simulation_time")) 
    {
        pconfig->final_time = strtod(value, NULL);
        pconfig->final_time_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "sigma")) 
    {
        pconfig->sigma = strtod(value, NULL);
        pconfig->sigma_was_set = true;
    }
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "beta")) 
    {
        pconfig->beta = strtod(value, NULL);
        pconfig->beta_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "cm")) 
    {
        pconfig->cm = strtod(value, NULL);
        pconfig->cm_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "print_rate")) 
    {
        pconfig->print_rate = (int)strtol (value, NULL, 10);
        pconfig->print_rate_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "output_dir")) 
    {
        pconfig->out_dir_name = strdup(value);
        pconfig->out_dir_name_was_set = true;
    }
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "output_steady_state_dir")) 
    {
        pconfig->out_steady_state_dir = strdup(value);
        pconfig->out_steady_state_dir_was_set = true;
    }
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "steady_state_filename")) 
    {
        pconfig->steady_state_filename = strdup(value);
        pconfig->steady_state_filename_was_set = true;
    }
    else if(SECTION_STARTS_WITH(STIM_SECTION)) 
    {


        if( pconfig->stim_configs == NULL) 
        {
            pconfig->stim_configs = new_stim_config();
        }

        struct stim_config *tmp = pconfig->stim_configs;

        if(MATCH_NAME("stim_start")) 
        {
            tmp->stim_start = (double)strtod(value, NULL);
            tmp->stim_start_was_set = true;
        } 
        else if(MATCH_NAME("stim_duration")) 
        {
            tmp->stim_duration = (double)strtod(value, NULL);
            tmp->stim_duration_was_set = true;
        }
        else if(MATCH_NAME("stim_current")) 
        {
            tmp->stim_current = (double)strtod(value, NULL);
            tmp->stim_current_was_set = true;
        }
        else if(MATCH_NAME("n_cycles")) 
        {
            tmp->n_cycles = (int)strtol(value, NULL, 10);
            tmp->n_cycles_was_set = true;
        }
        else if(MATCH_NAME("start_period")) 
        {
            tmp->start_period = (double)strtod(value, NULL);
            tmp->start_period_was_set = true;
        }
        else if(MATCH_NAME("end_period")) 
        {
            tmp->end_period = (double)strtod(value, NULL);
            tmp->end_period_was_set = true;
        }   
        else if(MATCH_NAME("period_step")) 
        {
            tmp->period_step = (double)strtod(value, NULL);
            tmp->period_step_was_set = true;
        }
        else if(MATCH_NAME("id_limit")) 
        {
            tmp->id_limit = (int)strtol(value, NULL, 10);
            tmp->id_limit_was_set = true;
        } 
        else 
        {
            //name is a reserved word in stim config
            if(MATCH_NAME("name")) 
            {
                fprintf(stderr, "name is a reserved word and should not be used inside a stimulus config section. Found in %s. Exiting!\n", section);
                exit(EXIT_FAILURE);
            }
        }
    }
    else if(MATCH_SECTION(PURKINJE_SECTION)) 
    {


        if(pconfig->purkinje_config == NULL) 
        {
            pconfig->purkinje_config = new_purkinje_config();
        }

        struct purkinje_config *tmp = pconfig->purkinje_config;

        if (MATCH_NAME ("start_discretization")) 
        {
            tmp->start_h = strtod (value, NULL);
            tmp->start_h_was_set = true;
        }
        else if(MATCH_NAME("name")) 
        {
            tmp->name = strdup(value);
            tmp->name_was_set = true;

        }
        else if(MATCH_NAME("network_file")) 
        {
            tmp->network_filename = strdup(value);
            tmp->network_filename_was_set = true;
        }
    }
    else 
    {
        fprintf(stderr, "Invalid name %s in section %s on the config file!\n", name, section);
        return 0;
    }

    return 1;
}