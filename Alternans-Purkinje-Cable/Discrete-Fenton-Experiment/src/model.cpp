#include "../include/model.h"

void Usage (const char pName[])
{
    printf("=====================================================================================\n");
    printf("Usage:> %s <input_filename>\n",pName);
    printf("-------------------------------------------------------------------------------------\n");
    printf("<input_filename> = Input filename with the parameters of the simulation\n");
    printf("=====================================================================================\n");
    /*
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
    */
}

void solveModel (int argc, char *argv[])
{
    Model *model = new Model(argc,argv);
    model->solve();
    //delete model;
}

Model::Model (int argc, char *argv[])
{
    User_Options *options = new User_Options(argc,argv);
    //options->print_user_options();

    if (options->steady_state)
    {
        sst = new SteadyState(options);
        sol = NULL;
    }
    else
    {
        sol = new Solver(options);
        sst = NULL;
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