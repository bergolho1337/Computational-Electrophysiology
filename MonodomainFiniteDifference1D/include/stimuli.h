#ifndef STIMULI_H
#define STIMULI_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

double get_spatial_stim_currents (const double x);
void compute_stimulus (double *stims, const double cur_time, const int np, const double dx);

#endif
