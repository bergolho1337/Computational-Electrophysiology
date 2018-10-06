#include <iostream>
#include "graph/graph.h"
#include "monodomain/ode_solver.h"
#include "monodomain/monodomain_solver.h"
#include "config/user_config.h"
#include "purkinje/purkinje.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        display_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    else
    {
        struct user_options *the_options = new_user_options(argc,argv);

        struct graph *the_purkinje_network = new_graph();

        struct monodomain_solver *the_monodomain_solver = new_monodomain_solver();

        struct ode_solver *the_ode_solver = new_ode_solver();

        configure_purkinje_network_from_options (the_purkinje_network,the_options);
        configure_monodomain_from_options (the_monodomain_solver,the_options);
        configure_ode_solver(the_ode_solver,the_purkinje_network->total_nodes);

        free_graph(the_purkinje_network);
        free_user_options(the_options);
    }
     
    return 0;
}