#include "noble.h"

/* ====================================================================================================== */
// dV/dt

double I_Stim__Nob (int k, double t)
{
  return 0;
}

double I_Leak__Nob (double t, double *y)
{
	return ((g_L*(y[0]-E_L)));
}

// Oscilacoes ainda presentes podem ser eliminadas somando um termo 1.0e+02
double g_K2__Nob (double t, double *y)
{
  return ((1.2e+03*pow(y[3],4.0e+00)));
	//return ((1.2e+03*pow(y[3],4.0e+00))+1.0e+02);          // Sem oscilacoes
}

// Trocar o primeiro coeficiente de 1.2e+03 para 1.3e+03 tambem tira o caracter auto-oscilatorio
double g_K1__Nob (double t, double *y)
{
	return (((1.2e+03*exp((((-y[0])-9.0e+01)/5.0e+01)))+(1.5e+01*exp(((y[0]+9.0e+01)/6.0e+01)))));
  //return (((1.3e+03*exp((((-y[0])-9.0e+01)/5.0e+01)))+(1.5e+01*exp(((y[0]+9.0e+01)/6.0e+01)))));    // Sem oscilacoes
}

double I_K__Nob (double t, double *y)
{
  return (((g_K1__Nob(t,y)+g_K2__Nob(t,y))*(y[0]+1.0e+02)));
}

double g_Na__Nob (double t, double *y)
{
  return ((pow(y[1],3.0e+00)*y[2]*g_Na_max));
}

// Mudando o coeficiente de 1.4e+02 para 1.225e+02 pode-se eliminar o comportamento auto-oscilatorio
double I_Na__Nob (double t, double *y)
{
	return ((g_Na__Nob(t,y)+1.4e+02)*(y[0]-E_Na));
  //return ((g_Na(t,y)+1.225e+02)*(y[0]-E_Na));
}

double dvdt__Nob (int k, double t, double *y)
{
	//return ((-(I_Na(t,y)+I_K(t,y)+I_Leak(t,y))/CM));
  return ((-(I_Na__Nob(t,y)+I_K__Nob(t,y)+I_Leak__Nob(t,y))+I_Stim__Nob(k,t))/CM);
}

/* ====================================================================================================== */
// dm/dt

double beta_m__Nob (double t, double *y)
{
	return (((1.2e+02*(y[0]+8.0e+00))/(exp(((y[0]+8.0e+00)/5.0e+00))-1.0e+00)));
}

double alpha_m__Nob (double t, double *y)
{
	return (((1.0e+02*((-y[0])-4.8e+01))/(exp((((-y[0])-4.8e+01)/1.5e+01))-1.0e+00)));
}

double dmdt__Nob (int k, double t, double *y)
{
	return ((alpha_m__Nob(t,y)*(1.0e+00-y[1]))-(beta_m__Nob(t,y)*y[1]));
}

/* ====================================================================================================== */
// dh/dt

double beta_h__Nob (double t, double *y)
{
	return ((1.0e+03/(1.0e+00+exp((((-y[0])-4.2e+01)/1.0e+01)))));
}

double alpha_h__Nob (double t, double *y)
{
	return ((1.7e+02*exp((((-y[0])-9.0e+01)/2.0e+01))));
}

double dhdt__Nob (int k, double t, double *y)
{
	return ((alpha_h__Nob(t,y)*(1.0e+00-y[2]))-(beta_h__Nob(t,y)*y[2]));
}

/* ====================================================================================================== */
// dn/dt

double beta_n__Nob (double t, double *y)
{
	return ((2.0e+00*exp((((-y[0])-9.0e+01)/8.0e+01))));
}

double alpha_n__Nob (double t, double *y)
{
	return (((1.0e-01*((-y[0])-5.0e+01))/(exp((((-y[0])-5.0e+01)/1.0e+01))-1.0e+00)));
}

double dndt__Nob (int k, double t, double *y)
{
	return ((alpha_n__Nob(t,y)*(1.0e+00-y[3]))-(beta_n__Nob(t,y)*y[3]));
}

/* ====================================================================================================== */

void setInitialConditions__Nob (double *y, int num_eq)
{
  y[0] = v0__Nob;
  y[1] = m0__Nob;
  y[2] = h0__Nob;
  y[3] = n0__Nob;
}

void setFunctions__Nob (Func *f, int num_eq)
{
  f[0] = dvdt__Nob;
  f[1] = dmdt__Nob;
  f[2] = dhdt__Nob;
  f[3] = dndt__Nob;
}
