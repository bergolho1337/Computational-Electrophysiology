#include "../include/sst.h"

// Number of threads to solve the system of ODEs
static constexpr int nthreads = 2;

SteadyState::SteadyState (int argc, char *argv[])
{
    dt = atof(argv[2]);
    tmax = atof(argv[3]);
    mesh_filename = argv[4];
    steady_filename = argv[5];
    M = nearbyint(tmax/dt);
    setSensibilityParam(argc,argv);
    g = new Graph(mesh_filename,dx);

    setControlVolumes();
    setFunctions();
    setInitCondModel();

    //g->printterm();
    //g->print();
    //print();
}

void SteadyState::solve ()
{
    FILE *sstFile = fopen(steady_filename.c_str(),"w+");
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

    #ifdef OUTPUT
    printf("[!] Solving transient problem ... \n");
    printf("[!] Progress\n");
    fflush(stdout);
    #endif

    // Time loop
    for (int i = 0; i < M; i++)
    {
        double t = i*dt;

        #ifdef OUTPUT
        printProgress(i,M);
        #endif

        //if (i % 200 == 0) writeVTKFile(i);
        if (i == 60000) writeSteadyStateFile(sstFile);

        // Solve the PDE 
        assembleLoadVector(b);
        x = sparseSolver.solve(b);
        moveVstar(x);

        // Solve the ODE
        solveODE(t);

        nextTimestep();
    }
    fclose(sstFile);

    #ifdef OUTPUT
    printf("ok\n");
    #endif
}

void SteadyState::moveVstar (const VectorXd vm)
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

void SteadyState::solveODE (double t)
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

void SteadyState::print ()
{
    printf("--------- Steady State -------------\n");
    printf("dt = %lf\n",dt);
    printf("tmax = %lf\n",tmax);
    printf("tmax/dt = %d\n",M);
    printf("dx = %lf\n",dx);
    printf("mesh_filename = %s\n",mesh_filename.c_str());
    printf("steady_filename = %s\n",steady_filename.c_str());
    printf("alfa = %lf\n",alfa);
    printf("d1 = %lf\n",d1);
    printf("------------------------------------\n");
}

// Build the coefficient matrix considering cells with the same diameter
void SteadyState::setMatrix (SpMat &a)
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
                    coeff.push_back(T(u,v,value));
                    ptrl = ptrl->next;
                }
                value = (ptr->num_edges*SIGMA*dx) + ALPHA;
                coeff.push_back(T(u,u,value));
            }
            // Is a special link to a Purkinje cell - PMJ
            else
            {
                double sum = C;
                while (ptrl != NULL)
                {
                    int v = ptrl->dest->id;
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

// Build the coefficient matrix considering cells with a different diameter
void SteadyState::setMatrix2 (SpMat &a)
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

void SteadyState::assembleLoadVector (VectorXd &b)
{
    double ALPHA = (BETA*Cm*dx*dx*dx) / dt;

    int np = b.size();
    for (int i = 0; i < np; i++)
        b(i) = vol[i].yOld[0] * ALPHA;
}

void SteadyState::setSensibilityParam (int argc, char *argv[])
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
}

void SteadyState::setControlVolumes ()
{
    // Capture the number of equations of the celullar model
    Node *ptr = g->getListNodes();
    int neq = num_eq;
    int np = g->getTotalNodes();
    vol = (Volume*)malloc(sizeof(Volume)*np);
    for (int i = 0; i < np; i++, ptr = ptr->next)
    {
        vol[i].type = ptr->type;
        vol[i].yOld = (double*)calloc(neq,sizeof(double));
        vol[i].yNew = (double*)calloc(neq,sizeof(double));
        vol[i].yStar = (double*)calloc(neq,sizeof(double));
    }
}

void SteadyState::setFunctions ()
{
  func = (Func*)malloc(sizeof(Func)*num_eq);  
  func[0] = dvdt__Nob;
  func[1] = dmdt__Nob;
  func[2] = dhdt__Nob;
  func[3] = dndt__Nob;
}

void SteadyState::setInitCondModel ()
{
    int neq = num_eq;
    int np = g->getTotalNodes();
    for (int i = 0; i < np; i++)
        for (int j = 0; j < neq; j++)
            vol[i].yOld[j] = y0__Nob[j];
}

bool SteadyState::isConnToPMJ (Edge *e)
{
    while (e != NULL)
    {
        if (e->dest->type == 1) return true;
        e = e->next;
    }
    return false;
}

void SteadyState::nextTimestep ()
{
    int np = g->getTotalNodes();
    for (int i = 0; i < np; i++) swap(&vol[i].yOld,&vol[i].yNew);
}

void swap (double **a, double **b)
{
    double *tmp = *a;
    *a = *b;
    *b = tmp;
}

void printProgress (int iter, int max_iter)
{
    double progress = iter / (double)max_iter;
    int barWidth = 100;

    cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) 
    {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    cout << "] " << int(progress * 100.0) << " %\r";
    cout.flush();
}

void SteadyState::writeVTKFile (int iter)
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

void SteadyState::writeSteadyStateFile (FILE *sstFile)
{
    int neq = num_eq;
    int np = g->getTotalNodes();
    for (int i = 0; i < np; i++)
    {
        for (int j = 0; j < neq; j++)
            fprintf(sstFile,"%.10lf ",vol[i].yOld[j]);
        fprintf(sstFile,"\n");
    }
}