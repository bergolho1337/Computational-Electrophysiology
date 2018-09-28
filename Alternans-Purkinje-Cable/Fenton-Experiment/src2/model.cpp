#include "../include/model.h"

void Usage (const char pName[])
{
  printf("=====================================================================================\n");
  printf("Usage:> %s -t/s <dt> <tmax> <mesh_file> <steady_state_file> <plot_file> [ALPHA] [d1] [SIGMA]\n",pName);
  printf("-------------------------------------------------------------------------------------\n");
  printf("-t = Steady-State mode\n");
  printf("-s = Solver mode\n\n");
  printf("<dt> = Size of the time discretization\n");
  printf("<tmax> = Maximum simulation time\n");
  printf("<mesh_file> = File with the mesh points and conections\n");
  printf("<steady_state_file> = File with the steady state solution\n");
  printf("<plot_file> = File with ids to be plotted\n");
  printf("-------------------------------------------------------------------------------------\n");
  printf(" !! Optional parameters !! (For the sensibility analysys)\n");
  printf("[ALPHA] = R_pmj*Vol_pmj\n");
  printf("[d1] = Diameter Purkinje cell\n");
  printf("[SIGMA] = Conductivity Gap Junction + Citoplasm\n");
  printf("-------------------------------------------------------------------------------------\n");
  printf("Steady-State Example: %s -t 0.1 5000 cable_dog.msh steady_dog.sst cable_dog.plt 1.375 0.002 0.004\n",pName);
  printf("Solver Example: %s -s 0.1 1000 cable_dog.msh steady_dog.sst cable_dog.plt 1.375 0.002 0.004\n",pName);
  printf("-------------------------------------------------------------------------------------\n");
}

void solveModel (int argc, char *argv[])
{
    Model *model = new Model(argc,argv);
    model->solve();
    //delete model;
}

Model::Model (int argc, char *argv[])
{
    string mode(argv[1]);
    if (mode == "-t")
    {
        sst = new SteadyState(argc,argv);
        sol = NULL;
    }
    else if (mode == "-s")
    {
        sol = new Solver(argc,argv);
        sst = NULL;
    }
    else
    {
        error("Invalid mode");
    }
}

void Model::solve ()
{
    if (sst != NULL)
        sst->solve();
    else
        sol->solve();
}

void Model::error (const char msg[])
{
    printf("[-] ERROR on Model ! %s !\n",msg);
    exit(EXIT_FAILURE);
}