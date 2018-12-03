#ifndef MONODOMAIN_H
#define MONODOMAIN_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "celular_model.h"

void solve_diffusion (double *sv, const double dx, const double dt,\
 			const int np, const int nodes,\
			double *vm);
void update_state_vector (double *sv, const double *vm,\
			  const int np, const int nodes);
void solve_reaction (double *sv, double *stims, const double t,\
		     const int np, const int nodes, const double dt);

#endif
