#include "../include/fitz.h"

void setInitialConditions__Fitz (double *y, int num_eq)
{
    y[0] = v0__Fitz;
    y[1] = w0__Fitz;
}

void setFunctions__Fitz (Func *f, int num_eq)
{
    f[0] = dvdt__Fitz;
    f[1] = dwdt__Fitz;
}

double I_stim__Fitz (double t)
{
    if (t >= 0 && t < 2)
        return v_stim__Fitz;
    else
        return 0;
}

double dvdt__Fitz (int k, double t, double *y)
{
    return y[0] - pow(y[0],3)/3.0 - y[1] + I_stim__Fitz(t);
}

double dwdt__Fitz (int k, double t, double *y)
{
    return (y[0] + a__Fitz - b__Fitz*y[1]) / tal__Fitz;
}