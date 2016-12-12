#include "hodkin.h"

// 0 = Vm,
// 1 = n,
// 2 = m,
// 3 = h

/* ====================================================================================================== */
// dV/dt

double I_Stim__Hod (int k, double t)
{
  if (t < 1)
    return v_stim__Hod;
  else
    return 0;
}

double I_Leak__Hod (double t, double *y)
{
  return g_L__Hod*(y[0]-E_L__Hod);
}

double I_K__Hod (double t, double *y)
{
  return g_K__Hod*pow(y[1],4)*(y[0]-E_K__Hod);
}

double I_Na__Hod (double t, double *y)
{
  return g_Na__Hod*pow(y[2],3)*y[3]*(y[0]-E_Na__Hod);
}

double dvdt__Hod (int k, double t, double *y)
{
  return (( I_Stim__Hod(k,t) - I_Na__Hod(t,y) - I_K__Hod(t,y) - I_Leak__Hod(t,y) ) / CM__Hod) ;
}

/* ====================================================================================================== */
// dn/dt

double beta_n__Hod (double t, double *y)
{
  return 0.125 * exp( (-y[0]-65.0) / 80.0);
}

double alpha_n__Hod (double t, double *y)
{
  return ( 0.01 * (y[0] + 55.0) ) / (1.0 - exp( (-y[0]-55.0) / 10.0 ) );
}

double dndt__Hod (int k, double t, double *y)
{
  return (alpha_n__Hod(t,y)*(1-y[1])) - (beta_n__Hod(t,y)*y[1]);
}

/* ====================================================================================================== */
// dm/dt

double beta_m__Hod (double t, double *y)
{
  return 4.0 * exp((-y[0]-65.0)/18.0);
}

double alpha_m__Hod (double t, double *y)
{
  return ( 0.1 * (y[0]+40.0) ) / ( 1.0 - exp( (-y[0]-40.0) / 10.0 ) );
}

double dmdt__Hod (int k, double t, double *y)
{
  return (alpha_m__Hod(t,y)*(1-y[2])) - (beta_m__Hod(t,y)*y[2]);
}

/* ********************************************************************************************************************** */
// dh/dt

double beta_h__Hod (double t, double *y)
{
  return 1.0 / ( exp((-y[0]-35.0)/10.0) + 1.0 );
}

double alpha_h__Hod (double t, double *y)
{
  return 0.07 * exp((-y[0]-65.0)/20.0 );
}

double dhdt__Hod (int k, double t, double *y)
{
  return (alpha_h__Hod(t,y)*(1-y[3])) - (beta_h__Hod(t,y)*y[3]);
}

/* ====================================================================================================== */

void setInitialConditions__Hod (double *y, int num_eq)
{
  memset(y,0,sizeof(double)*num_eq);
  y[1] = alpha_n__Hod(0,y)/(alpha_n__Hod(0,y)+beta_n__Hod(0,y));
  y[2] = alpha_m__Hod(0,y)/(alpha_m__Hod(0,y)+beta_m__Hod(0,y));
  y[3] = alpha_h__Hod(0,y)/(alpha_h__Hod(0,y)+beta_h__Hod(0,y));
  //y[0] = v0__Hod;
  //y[1] = n0__Hod;
  //y[2] = m0__Hod;
  //y[3] = h0__Hod;
  printf("v0 = %e\n",y[0]);
  printf("n0 = %e\n",y[1]);
  printf("m0 = %e\n",y[2]);
  printf("h0 = %e\n",y[3]);
}

void setFunctions__Hod (Func *f, int num_eq)
{
  f[0] = dvdt__Hod;
  f[1] = dndt__Hod;
  f[2] = dmdt__Hod;
  f[3] = dhdt__Hod;
}
