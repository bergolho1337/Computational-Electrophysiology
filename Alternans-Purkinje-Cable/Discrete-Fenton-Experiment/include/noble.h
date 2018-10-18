#ifndef NOBLE_H_
#define NOBLE_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

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
static const int num_eq = 4;

// Stimulus configuration
static const double stim_current = 1000.0f;
static const double stim_start = 0.0f;
static const double stim_duration = 2.0f;
static const double start_period = 300.0f;
static const double end_period = 300.0f;
static const double period_step = 100.0f;
static const int n_cycles = 20;
static const int id_limit = 20;

// Cycle length
static const double cycle_length = 300.0;  

// Initial conditions (default)
static const double y0__Nob[4] = {-75.5344986658,0.0605467272,0.7259001355,0.4709239708};

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

#endif