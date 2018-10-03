#include "../include/solver.h"

// Number of threads to solve the system of ODEs
static constexpr int nthreads = 2;

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
    //setPlot(); 
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
    int max_num_digits = get_num_digits(M);

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
        //writePlotData(t);
        if (i % 10 == 0)
            writeStateVectorToFile(i,max_num_digits,t);

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
        d1 = 0.003;
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
    //BETA = 0.14;
    //cout << "Beta = " << BETA << endl;
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
    vel->velocityFile = fopen("output/velocity.txt","w+");
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
        sprintf(filename,"output/data%d.dat",plot->ids[i]);
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

    double ALPHA = (BETA*Cm*dx*dx*dx) / dt;

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
                double value = -SIGMA * dx;
                while (ptrl != NULL)
                {
                    int v = ptrl->dest->id;
                    //printf("(%d,%d) = %.10lf\n",u,v,value);
                    coeff.push_back(T(u,v,value));
                    ptrl = ptrl->next;
                }
                value = (ptr->num_edges*SIGMA*dx) + ALPHA;
                //printf("(%d,%d) = %.10lf\n",u,u,value);
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
    
    // Print non-zero coefficients
    FILE *file = fopen("matrix.txt","w+");
    for (int i = 0; i < coeff.size(); i++)
        fprintf(file,"(%d,%d) = %.10lf\n",coeff[i].row(),coeff[i].col(),coeff[i].value());
    fclose(file);

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
    double ALPHA = (BETA*Cm*dx*dx*dx) / dt;

    int np = b.size();
    for (int i = 0; i < np; i++)
        b(i) = vol[i].yOld[0] * ALPHA;
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

void Solver::writeStateVectorToFile (int count, int max_num_digits, double t)
{
    int num_digits = get_num_digits(count);
    int num_zeros = max_num_digits - num_digits;

    string str_zeros = "";
    for (int i = 0; i < num_zeros; i++)
        str_zeros += "0";
    
    stringstream ss;
    ss << "tmp/sv_" << str_zeros << count;
    FILE *file = fopen(ss.str().c_str(),"w+");
    int neq = num_eq;
    int ncells = g->getTotalNodes();

    for (int i = 0; i < ncells; i++)
    {
        fprintf(file,"%.10lf ",t);
        for (int j = 0; j < num_eq-1; j++)
        {
            fprintf(file,"%.10lf ",vol[i].yOld[j]);
        }
        fprintf(file,"%.10lf\n",vol[i].yOld[(num_eq-1)]);
    }

    fclose(file);
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
    FILE *vFile = fopen("output/v.txt","w+");
    FILE *dFile = fopen("output/delay.txt","w+");
    FILE *tFile = fopen("output/time.txt","w+");
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

int get_num_digits (int num)
{
    if (num == 0) return 1;

    int num_digits = 0;
    while (num != 0)
    {
        num /= 10;
        num_digits++;
    }
    return num_digits;
}