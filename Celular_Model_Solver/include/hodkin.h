#ifndef HODKIN_H_
#define HODKIN_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// Function pointer
typedef double (*Func) (int k, double t, double *y);

/* ============================== CONSTANTS ========================================== */
// Initial conditions
const double v0__Hod = 0.0;
//const double n0__Hod = 0.32;
//const double m0__Hod = 0.05;
//const double h0__Hod = 0.6;
// Parameters:
const double E_Na__Hod = 50.0;             // mV
const double g_Na__Hod = 120.0;              // mS/cm^2
const double E_K__Hod = -77.0;                // mV
const double g_K__Hod = 36.0;                // mS/cm^2
const double E_L__Hod = -54.387;                // mV
const double g_L__Hod = 0.3;                 // mS/cm^2
const double CM__Hod = 1.0;                  // uF/cm^2
const double v_stim__Hod = 0.0;              // uA/cm^2

/* ********************************************************************************************************************** */
// dV/dt

double I_Stim__Hod (int k, double t);
double I_Leak__Hod (double t, double *y);
double I_K__Hod (double t, double *y);
double I_Na__Hod (double t, double *y);
double dvdt__Hod (int k, double t, double *y);
/* ********************************************************************************************************************** */
// dm/dt

double beta_m__Hod (double t, double *y);
double alpha_m__Hod (double t, double *y);
double dmdt__Hod (int k, double t, double *y);
/* ********************************************************************************************************************** */
// dh/dt

double beta_h__Hod (double t, double *y);
double alpha_h__Hod (double t, double *y);
double dhdt__Hod (int k, double t, double *y);
/* ********************************************************************************************************************** */
// dn/dt

double beta_n__Hod (double t, double *y);
double alpha_n__Hod (double t, double *y);
double dndt__Hod (int k, double t, double *y);
/* ********************************************************************************************************************** */

void setInitialConditions__Hod (double *y, int num_eq);
void setFunctions__Hod (Func *f, int num_eq);

#endif
