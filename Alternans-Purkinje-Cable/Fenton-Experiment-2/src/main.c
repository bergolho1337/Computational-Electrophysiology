#include "alg/grid/grid.h"
#include "ini_parser/ini.h"
#include "monodomain/monodomain_solver.h"
#include "monodomain/ode_solver.h"
#include "utils/logfile_utils.h"
#include "monodomain/output_utils.h"
#include "string/sds.h"

int main (int argc, char **argv) 
{

    struct user_options *options;
    options = new_user_options ();

    struct grid *the_grid;
    the_grid = new_grid ();

    struct monodomain_solver *monodomain_solver;
    monodomain_solver = new_monodomain_solver ();

    struct ode_solver *ode_solver;
    ode_solver = new_ode_solver ();

    // ______________________________________________________________________________________
    // Configuration step
    // First we have to get the config file path
    get_config_file (argc, argv, options);

    if (options->config_file) 
    {
        // Here we parse the config file
        if (ini_parse (options->config_file, parse_config_file, options) < 0) 
        {
            fprintf (stderr, "Error: Can't load the config file %s\n", options->config_file);
            return EXIT_FAILURE;
        }
    }

    // The command line options always overwrite the config file
    parse_options (argc, argv, options);

    // Create the output dir and the logfile
    if (options->out_dir_name) 
    {
        sds buffer = sdsnew ("");
        create_dir_if_no_exists (options->out_dir_name);
        buffer = sdscatfmt (buffer, "%s/outputlog.txt", options->out_dir_name);
        open_logfile (buffer);
        sdsfree (buffer);
    }
    if (options->out_steady_state_dir) 
    {
        create_dir_if_no_exists (options->out_steady_state_dir);
    }

    configure_ode_solver_from_options (ode_solver, options);
    configure_monodomain_solver_from_options (monodomain_solver, options);
    configure_grid_from_options (the_grid, options);


    #ifndef COMPILE_CUDA
    if (ode_solver->gpu) 
    {
        print_to_stdout_and_file ("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
        ode_solver->gpu = false;
    }
    #endif

    int np = monodomain_solver->num_threads;

    if (np == 0)
        np = 1;

    #if defined(_OPENMP)
    omp_set_num_threads (np);
    #endif

    // ______________________________________________________________________________________
    // Solving the problem
    solve_monodomain (monodomain_solver, ode_solver, the_grid, options);

    clean_and_free_grid (the_grid);
    free_ode_solver (ode_solver);

    free (monodomain_solver);

    free_user_options (options);
    close_logfile ();

    return EXIT_SUCCESS;
}
