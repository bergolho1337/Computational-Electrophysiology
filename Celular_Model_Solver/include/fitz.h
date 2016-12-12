#ifndef FITZ_H_
#define FITZ_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// Maximum time = 200s
// dt = 0.1

// Function pointer
typedef double (*Func) (int k, double t, double *y);

/* ============================== CONSTANTS ========================================== */
const double v0__Fitz = 0.0;
const double w0__Fitz = 0.0;
const double v_stim__Fitz = 0.1;
const double a__Fitz = 0.0;                         // Auto-oscillatory
const double b__Fitz = 0.0;                         // Auto-oscillatory
const double tal__Fitz = 12.5;
/* ============================ FUNCTIONS ============================================= */
double dvdt__Fitz (int k, double t, double *y);
double dwdt__Fitz (int k, double t, double *y);
double I_stim__Fitz (double t);
/* ==================================================================================== */

void setInitialConditions__Fitz (double *y, int num_eq);
void setFunctions__Fitz (Func *f, int num_eq);

#endif
