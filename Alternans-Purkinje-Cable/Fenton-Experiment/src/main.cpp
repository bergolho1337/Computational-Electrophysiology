/*
========= FINITE VOLUME METHOD ----- MONODOMAIN EQUATION (CABLE) =============================================
  Problem: solve the cable equation using the monodomain equation applied to the FVM
    { BETA*Cm*V_t = SIGMA*V_xx
    { V'(0,t) = V'(1,t) = I_PMJ           Non-homogeneus Neumann condition at the terminal volumes
    { V(x,0) = V_inf
***********************************************************************************************************
    BETA = Ratio surface-area per volume of the cell (cm^-1)
    Cm = Membrane capacitance of the celullar membrane (uF/cm^2)
    SIGMA = Conductivity of the celullar membrane (mS/cm^2) -- All the fiber has the same conductivity
***********************************************************************************************************
==============================================================================================================
*/

#include <iostream>
#include <string>
#include "monodomain/monodomain_solver.h"
#include "monodomain/ode_solver.h"
#include "ini_parser/ini.h"
#include "monodomain/output_utils.h"
#include "utils/logfile_utils.h"

using namespace std;

int main (int argc, char *argv[])
{
  struct user_options *options;
  options = new_user_options ();

  struct grid *grid;
  grid = new_grid ();

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
      create_dir_if_no_exists (options->out_dir_name);
  }
  if (options->out_steady_state_dir) 
  {
      create_dir_if_no_exists (options->out_steady_state_dir);
  }

  configure_ode_solver_from_options (ode_solver, options);
  configure_monodomain_solver_from_options (monodomain_solver, options);
  configure_grid_from_options (grid, options);

  #ifndef COMPILE_CUDA
  if (ode_solver->gpu) 
  {
      print_to_stdout_and_file ("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
      ode_solver->gpu = false;
  }
  #endif

  // ______________________________________________________________________________________
  // Solving the problem
  solve_monodomain (monodomain_solver, ode_solver, grid, options);

  cout << "Hello world" << endl;
  return 0;
}
