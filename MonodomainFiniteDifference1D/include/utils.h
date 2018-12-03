#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#define PRINT_LINE "===================================================================================================="

void write_VTK_to_file (const double *sv, const double dx,\
                        const int Ncell, const int Nodes,\
			const int iter);
void write_plot_data (FILE *file, const double t, const double *sv,\
			const int Ncell, const int Nodes,\
			const int id);
void print_stimulus (const double *stim_current, const int Ncell, const double dx);
void print_state_vector (const double *sv, const int Ncell, const int Nodes);
void print_progress (int iter, int max_iter);
void print_configuration_parameters(const double dx, const double dt, const double tmax, const double lmax,\
					const int Ncell, const int Niter, const int Nodes,\
					const int plot_cell_id, const int print_rate);


#endif
