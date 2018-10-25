#include "../include/model.h"

void Usage (const char pName[])
{
    printf("=====================================================================================\n");
    printf("Usage:> %s <input_filename>\n",pName);
    printf("-------------------------------------------------------------------------------------\n");
    printf("<input_filename> = Input filename with the parameters of the simulation\n");
    printf("=====================================================================================\n");
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