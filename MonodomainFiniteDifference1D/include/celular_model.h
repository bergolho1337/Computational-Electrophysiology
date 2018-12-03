#ifndef CELULAR_MODEL_H
#define CELULAR_MODEL_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

void compute_initial_conditions (double *sv, const int Ncell, const int Nodes);
double dvdt (const double V, const double m, const double h, const double n, const double stim_current);
double dmdt (const double V, const double m);
double dhdt (const double V, const double h);
double dndt (const double V, const double n);

#endif
