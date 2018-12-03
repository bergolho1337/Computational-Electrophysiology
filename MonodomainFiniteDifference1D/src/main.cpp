// --------------------------------------------------------------------------------------------------------------------
// This program solves the 1D Monodomain equation using the Mitchell-Shaeffer celular model
// --------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>

#include "../include/timer.h"
#include "../include/utils.h"
#include "../include/stimuli.h"
#include "../include/monodomain.h"


using namespace std;

int main (int argc, char *argv[])
{
	// Print configuration
	FILE *file = fopen("output/sv.dat","w+");
	const int plot_cell_id = 250;	
	const int print_rate = 100;

	// Finite difference parameters configuration
	const double dx = 0.01;
	const double dt = 0.01;
	const double tmax = 3000.0;
	const double lmax = 5.0;

	int Ncell = nearbyint(lmax / dx);
	int Niter = nearbyint(tmax / dt);

	cout << "Number of cells = " << Ncell << endl;
	cout << "Number of iterations = " << Niter << endl;
	
	// OpenMP configuration
	omp_set_dynamic(0);
    	omp_set_num_threads(1);

	// ODE's parameters
	int Nodes = 4;		// Noble, 1962

	print_configuration_parameters(dx,dt,tmax,lmax,\
					Ncell,Niter,Nodes,\
					plot_cell_id,print_rate);

	// Allocate memory
	double *sv = new double[Ncell*Nodes]();
	double *stim_current = new double[Ncell]();
	double *vm = new double[Ncell]();
	
	compute_initial_conditions(sv,Ncell,Nodes);
	//print_state_vector(sv,Ncell,Nodes);

	double start, finish, elapsed;
	GET_TIME(start);	

	for (int k = 0; k < Niter; k++)
	{
		double t = dt*k;

		print_progress(k,Niter);	

		if (k % print_rate == 0)
		{
			//write_VTK_to_file(sv,dx,Ncell,Nodes,k);
			write_plot_data(file,t,sv,Ncell, Nodes, plot_cell_id);
		}

		compute_stimulus(stim_current,t,Ncell,dx);
		//print_stimulus(stim_current,Ncell,dx);
		
		solve_diffusion(sv,dx,dt,Ncell,Nodes,vm);

		update_state_vector(sv,vm,Ncell,Nodes);		
		
		solve_reaction(sv,stim_current,t,Ncell,Nodes,dt);

	}
	GET_TIME(finish);
	elapsed = finish - start;
	printf("%s\n",PRINT_LINE);
	printf("Elapsed time = %.10lf\n",elapsed);
	printf("%s\n",PRINT_LINE);

	// Free memory
	delete [] sv;
	delete [] stim_current;
	delete [] vm;
	fclose(file);	

	return 0;  
}
