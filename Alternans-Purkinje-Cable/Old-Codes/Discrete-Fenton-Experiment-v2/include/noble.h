#ifndef NOBLE_H_
#define NOBLE_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "stimulus.h"

// Maximum time = 2000 ms
// dt = 0.1 ms
// y[0] = vm
// y[1] = m
// y[2] = h
// y[3] = n

// Function pointer
typedef double (*Func) (int type, int point, double t, double y[]);

/* ============================== CONSTANTS ========================================== */

// Number of equations
static const int NUM_EQ = 4;

// Initial conditions (default)
static const double YO__NOBLE[4] = {-75.5344986658,0.0605467272,0.7259001355,0.4709239708};

// Steady state
// -79.0666103401 0.0499158420 0.8049140578 0.2603110528
// -75.5344986658 0.0605467272 0.7259001355 0.4709239708

// Parameters:
const double g_Na_max = 4.0e+02;
const double E_Na = 4.0e+01;
const double g_L = 7.5e-02;
const double E_L = -6.0e+01;
const double CM = 1.2e+01;

/* ********************************************************************************************************************** */
// dV/dt

double I_Stim__Nob (int point, double t);
double I_Leak__Nob (double t, double vm, double m, double h, double n);
double g_K2__Nob (double t, double vm, double m, double h, double n);
double g_K1__Nob (double t, double vm, double m, double h, double n);
double I_K__Nob (double t, double vm, double m, double h, double n);
double g_Na__Nob (double t, double vm, double m, double h, double n);
double I_Na__Nob (double t, double vm, double m, double h, double n);
double I_Na_NoOscilation__Nob (double t, double vm, double m, double h, double n);
double dvdt__Nob (int type, int point, double t, double y[]);
/* ********************************************************************************************************************** */
// dm/dt

double beta_m__Nob (double t, double vm, double m, double h, double n);
double alpha_m__Nob (double t, double vm, double m, double h, double n);
double dmdt__Nob (int type, int point, double t, double y[]);
/* ********************************************************************************************************************** */
// dh/dt

double beta_h__Nob (double t, double vm, double m, double h, double n);
double alpha_h__Nob (double t, double vm, double m, double h, double n);
double dhdt__Nob (int type, int point, double t, double y[]);
/* ********************************************************************************************************************** */
// dn/dt

double beta_n__Nob (double t, double vm, double m, double h, double n);
double alpha_n__Nob (double t, double vm, double m, double h, double n);
double dndt__Nob (int type, int point, double t, double y[]);
/* ********************************************************************************************************************** */

void set_celular_model (Stimulus *stim_config);

#endif