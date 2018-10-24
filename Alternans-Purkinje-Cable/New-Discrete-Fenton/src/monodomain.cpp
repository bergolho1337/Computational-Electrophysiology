#include "../include/monodomain.h"

void Usage (const char pName[])
{
    printf("=====================================================================================\n");
    printf("Usage:> %s <input_filename>\n",pName);
    printf("-------------------------------------------------------------------------------------\n");
    printf("<input_filename> = Input filename with the parameters of the simulation\n");
    printf("=====================================================================================\n");
}

void solve_monodomain (int argc, char *argv[])
{
    Monodomain *monodomain = new Monodomain(argc,argv);

    monodomain->solve();
}

Monodomain::Monodomain (int argc, char *argv[])
{
    user_options = new Options(argc,argv);
    //user_options->print_user_options();

    solver = new Solver(user_options);
    //solver->print();

}

void Monodomain::solve ()
{
    cout << "[Solver] Solving monodomain equation" << endl;
    
    solver->solve();
}