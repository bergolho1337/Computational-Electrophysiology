#ifndef NOBLE_H_
#define NOBLE_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// Maximum time = 2s
// dt = 0.0001

// Function pointer
typedef double (*Func) (int k, double t, double *y);

/* ============================== CONSTANTS ========================================== */
// Initial conditions
const double v0__Nob = -8.0e+01;
const double m0__Nob = 1.0e-02;
const double h0__Nob = 8.0e-01;
const double n0__Nob = 1.0e-02;
// Parameters:
const double g_Na_max = 4.0e+05;
const double E_Na = 4.0e+01;
const double g_L = 7.5e+01;
const double E_L = -6.0e+01;
const double CM = 1.2e+01;

/* ********************************************************************************************************************** */
// dV/dt

double I_Stim__Nob (int k, double t);
double I_Leak__Nob (double t, double *y);
double g_K2__Nob (double t, double *y);
double g_K1__Nob (double t, double *y);
double I_K__Nob (double t, double *y);
double g_Na__Nob (double t, double *y);
// Changing the 1.4e+02 to 1.225e+02 we can eliminate the auto-oscilatory behaviour of the Noble model
double I_Na__Nob (double t, double *y);
double dvdt__Nob (int k, double t, double *y);
/* ********************************************************************************************************************** */
// dm/dt

double beta_m__Nob (double t, double *y);
double alpha_m__Nob (double t, double *y);
double dmdt__Nob (int k, double t, double *y);
/* ********************************************************************************************************************** */
// dh/dt

double beta_h__Nob (double t, double *y);
double alpha_h__Nob (double t, double *y);
double dhdt__Nob (int k, double t, double *y);
/* ********************************************************************************************************************** */
// dn/dt

double beta_n__Nob (double t, double *y);
double alpha_n__Nob (double t, double *y);
double dndt__Nob (int k, double t, double *y);
/* ********************************************************************************************************************** */

void setInitialConditions__Nob (double *y, int num_eq);
void setFunctions__Nob (Func *f, int num_eq);

#endif
