#include "../include/sst.h"

// Number of threads to solve the system of ODEs
static constexpr int nthreads = 2;

SteadyState::SteadyState (User_Options *options)
{
    dt = options->dt;
    tmax = options->tmax;
    mesh_filename = options->mesh_filename;
    steady_filename = options->sst_filename;
    M = nearbyint(tmax/dt);
    setSensibilityParam(options);
    
    dx = options->start_h / options->num_div_cell;
    g = setPurkinjeMeshFromFile(mesh_filename,dx);
    g->setGapJunctions(options->num_div_cell);
    
    setControlVolumes();
    setFunctions();
    setInitCondModel();

    //g->print();
}

void SteadyState::solve ()
{
    FILE *sstFile = fopen(steady_filename.c_str(),"w+");
    int np = g->get_total_nodes();

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

    double ALPHA = (BETA*Cm*dx*dx*dx) / dt;

    // Time loop
    for (int i = 0; i < M; i++)
    {
        double t = i*dt;

        #ifdef OUTPUT
        printProgress(i,M);
        #endif

        //if (i % 100 == 0) writeVTKFile(i);
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
    int ncells = g->get_total_nodes();
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
    double A = (4.0*GGAP*dx) / (M_PI*d1*d1);
    double B = (SIGMA);
    double C = (BETA*Cm*dx*dx) / (dt);

    // Non-zero coefficients
    vector<T> coeff;

    double diagonal_value;
    Node *ptr = g->get_list_nodes();
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
    
    a.setFromTriplets(coeff.begin(),coeff.end());
    a.makeCompressed();
}

void SteadyState::assembleLoadVector (VectorXd &b)
{
    double C = (BETA*Cm*dx*dx) / (dt);

    int np = b.size();
    for (int i = 0; i < np; i++)
        b(i) = vol[i].yOld[0] * C;
    
}

void SteadyState::setSensibilityParam (User_Options *options)
{
    alfa = options->alfa;
    d1 = options->diameter;
    SIGMA = options->sigma_c;
    GGAP = options->G_gap;

    BETA = 4.0 / d1 * 1.0e-04;
}

void SteadyState::setControlVolumes ()
{
    // Capture the number of equations of the celullar model
    Node *ptr = g->get_list_nodes();
    int neq = num_eq;
    int np = g->get_total_nodes();
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
    int np = g->get_total_nodes();
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
    int np = g->get_total_nodes();
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
    Node *ptr = g->get_list_nodes();
    np = g->get_total_nodes();
    ne = g->get_total_edges();

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
    ptr = g->get_list_nodes();
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
    ptr = g->get_list_nodes();
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
    int np = g->get_total_nodes();
    for (int i = 0; i < np; i++)
    {
        for (int j = 0; j < neq; j++)
            fprintf(sstFile,"%.10lf ",vol[i].yOld[j]);
        fprintf(sstFile,"\n");
    }
}