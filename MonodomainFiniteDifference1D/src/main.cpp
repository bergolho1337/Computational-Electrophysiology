// --------------------------------------------------------------------------------------------------------------------
// This program solves the 1D Monodomain equation using the Mitchell-Shaeffer celular model
// --------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <omp.h>

#include "../include/timer.h"
#include "../include/config.h"
#include "../include/utils.h"
#include "../include/plot.h"
#include "../include/stimuli.h"
#include "../include/monodomain.h"

using namespace std;

int main (int argc, char *argv[])
{
	if (argc-1 != 1)
	{
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	// Initialize all the structures 
	struct user_options *options;
	options = new_user_options();
	
	struct monodomain_solver *solver;
	solver = new_monodomain_solver();

	struct plot_config *plotter;
	plotter = new_plot_config();

	struct stim_config *stim;
	stim = new_stim_config();

	// Reading input configuration file
	read_input_file(options,argv[1]);

	// Parse the user input into the structures
	configure_solver_from_options(solver,options);
	configure_plot_from_options(plotter,solver,options);
	configure_stimulus_from_options(stim,options);

	// OpenMP configuration
	omp_set_dynamic(0);
    	omp_set_num_threads(1);

	print_configuration_parameters(solver->dx,solver->dt,solver->tmax,solver->lmax,\
					solver->Ncell,solver->Niter,Nodes,\
					plotter->plot_cell_ids,plotter->print_rate,plotter->sst_rate);

	// Call the solver function for the monodomain equation
	solve_monodomain(solver,stim,plotter);
	
	// Free memory
	free_user_options(options);
	free_monodomain_solver(solver);
	free_plot_config(plotter);
	free_stim_config(stim);

	return 0;  
}
