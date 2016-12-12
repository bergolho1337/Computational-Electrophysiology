#include "mitchell.h"

double dvdt__Mit (int k, double t, double *y)
{
  return (I_in__Mit(y) + I_out__Mit(y) + I_stim__Mit(t));
}

double I_in__Mit (double *y)
{
  return (y[1] * Cv(y) / tal_in);
}

double I_out__Mit (double *y)
{
  return (-y[0] / tal_out);
}

// >>>> Pacing can be configured here <<<<
double I_stim__Mit (double t)
{
  if (t < 1 || (t > 500.0 && t < 501.0))
    return v_stim;
  else
    return 0;
}

double Cv (double *y)
{
  return (pow(y[0],2)*(1-y[0]));
}

double dhdt__Mit (int k, double t, double *y)
{
  if (y[0] < v_gate)
    return ((1-y[1]) / tal_open);
  else
    return (-y[1] / tal_close);
}

void setInitialConditions__Mit (double *y, int num_eq)
{
  y[0] = v0__Mit;
  y[1] = h0__Mit;
}

void setFunctions__Mit (Func *f, int num_eq)
{
  f[0] = dvdt__Mit;
  f[1] = dhdt__Mit;
}
