#include <iostream>
#include "monodomain/config/config_parser.h"
#include "ini_parser/ini.h"
#include "grid/grid.h"
#include "monodomain/monodomain_solver.h"
#include "monodomain/timer.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 != 2)
    {
        display_usage(argv);
        return 1;
    }
    else
    {
        // Build the user configuration structure
        struct user_options *options;
        options = new_user_options ();
        
        struct grid *the_grid;
        the_grid = new_grid ();

        struct monodomain_solver *monodomain_solver;
        monodomain_solver = new_monodomain_solver ();

        // First we have to get the config file path
        get_config_file (argc, argv, options);

        if (options->config_file.size()) 
        {
            // Here we parse the config file
            if (ini_parse (options->config_file.c_str(), parse_config_file, options) < 0) 
            {
                fprintf (stderr, "Error: Can't load the config file %s\n", options->config_file.c_str());
                return EXIT_FAILURE;
            }
        }
        configure_monodomain_solver_from_options(monodomain_solver,options);
        configure_grid_from_options(the_grid,options);

        double start, finish, elapsed;
        GET_TIME(start);
        
        solve_monodomain(monodomain_solver,the_grid,options);
        
        GET_TIME(finish);
        elapsed = finish - start;
        printf("==========================================================\n");
        printf("[!] Time elapsed = %.10lf seconds\n",elapsed);
        printf("==========================================================\n");

        clean_and_free_grid(the_grid);
        free_monodomain_solver(monodomain_solver);
        free_user_options(options);        

        return 0;
    }
}