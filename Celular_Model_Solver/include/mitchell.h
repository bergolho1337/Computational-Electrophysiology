#ifndef MITCHELL_H_
#define MITCHELL_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// Maximum time = 500ms
// dt = 0.1

// Function pointer
typedef double (*Func) (int k, double t, double *y);

/* ============================== CONSTANTS ========================================== */
const double tal_in = 0.3;
const double tal_out = 6.0;
const double tal_open = 120.0;
const double tal_close = 150.0;
const double v_gate = 0.13;
const double v_stim = 0.1;
const double v0__Mit = 0.0;
const double h0__Mit = 1.0;
/* ============================ FUNCTIONS ============================================= */
double dvdt__Mit (int k, double t, double *y);
double I_in__Mit (double *y);
double I_out__Mit (double *y);
double I_stim__Mit (double t);
double Cv (double *y);

double dhdt__Mit (int k, double t, double *y);
/* ==================================================================================== */

void setInitialConditions__Mit (double *y, int num_eq);
void setFunctions__Mit (Func *f, int num_eq);

#endif
