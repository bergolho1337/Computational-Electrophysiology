#include "../include/cell.h"

int __nbeats;

Cell::Cell (int neq)
{
  yOld.assign(neq,0);
  yNew.assign(neq,0);
}

void Cell::setCell ()
{
  yOld[0] = v0;
  yOld[1] = m0;
  yOld[2] = h0;
  yOld[3] = n0;
}

// Euler explicito
void Cell::solve (int id, double t, double dt)
{
    yNew[0] = yOld[0] + dvdt(id,t,yOld.data())*dt;
    yNew[1] = yOld[1] + dmdt(id,t,yOld.data())*dt;
    yNew[2] = yOld[2] + dhdt(id,t,yOld.data())*dt;
    yNew[3] = yOld[3] + dndt(id,t,yOld.data())*dt;    
}

void Cell::swap ()
{
    yOld.swap(yNew);
}

void Cell::write (double t, FILE *out)
{
    fprintf(out,"%lf ",t);
    for (int i = 0; i < (int)yOld.size()-1; i++)
        fprintf(out,"%lf ",yOld[i]);
    fprintf(out,"%lf\n",yOld.back());
}

/* ====================================================================================================== */
// dV/dt
double Cell::I_Stim (int id, double t)
{
  if (id < 2)
  {
    for (int beat = 0; beat < __nbeats; beat++)
    {
      double taux = beat*PACING;
      if (t >= taux && t <= taux + T_STIM)
        return V_STIM;
    }
  }
  return 0;
}

double Cell::I_Leak (double t, double *y)
{
	return ((G_L*(y[0]-E_L)));
}

// Oscilacoes ainda presentes podem ser eliminadas somando um termo 1.0e+02
double Cell::g_K2 (double t, double *y)
{
  return ((1.2*pow(y[3],4.0e+00)));
	//return ((1.2e+03*pow(y[3],4.0e+00))+1.0e-01);          // Sem oscilacoes
}

// Trocar o primeiro coeficiente de 1.2e+03 para 1.3e+03 tambem tira o caracter auto-oscilatorio
double Cell::g_K1 (double t, double *y)
{
	return (((1.2*exp((((-y[0])-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((y[0]+9.0e+01)/6.0e+01)))));
  //return (((1.3*exp((((-y[0])-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((y[0]+9.0e+01)/6.0e+01)))));    // Sem oscilacoes
}

double Cell::I_K (double t, double *y)
{
  return (((g_K1(t,y)+g_K2(t,y))*(y[0]+1.0e+02)));
}

double Cell::g_Na (double t, double *y)
{
  return ((pow(y[1],3.0e+00)*y[2]*G_NA_MAX));
}

// Mudando o coeficiente de 1.4e+02 para 1.225e+02 pode-se eliminar o comportamento auto-oscilatorio
double Cell::I_Na (double t, double *y)
{
	return ((g_Na(t,y)+1.4e-01)*(y[0]-E_NA));
  //return ((g_Na(t,y)+1.220e-01)*(y[0]-E_Na));      // Sem oscilacoes
}

double Cell::dvdt (int k, double t, double *y)
{
	//return ((-(I_Na(t,y)+I_K(t,y)+I_Leak(t,y))/CM));
  return ((-(I_Na(t,y)+I_K(t,y)+I_Leak(t,y))+I_Stim(k,t))/CM);
}

/* ====================================================================================================== */
// dm/dt

double Cell::beta_m (double t, double *y)
{
	return (((1.2e-01*(y[0]+8.0e+00))/(exp(((y[0]+8.0e+00)/5.0e+00))-1.0e+00)));
}

double Cell::alpha_m (double t, double *y)
{
	return (((1.0e-01*((-y[0])-4.8e+01))/(exp((((-y[0])-4.8e+01)/1.5e+01))-1.0e+00)));
}

double Cell::dmdt (int k, double t, double *y)
{
	return ((alpha_m(t,y)*(1.0e+00-y[1]))-(beta_m(t,y)*y[1]));
}

/* ====================================================================================================== */
// dh/dt

double Cell::beta_h (double t, double *y)
{
	return ((1.0/(1.0e+00+exp((((-y[0])-4.2e+01)/1.0e+01)))));
}

double Cell::alpha_h (double t, double *y)
{
	return ((1.7e-01*exp((((-y[0])-9.0e+01)/2.0e+01))));
}

double Cell::dhdt (int k, double t, double *y)
{
	return ((alpha_h(t,y)*(1.0e+00-y[2]))-(beta_h(t,y)*y[2]));
}

/* ====================================================================================================== */
// dn/dt

double Cell::beta_n (double t, double *y)
{
	return ((2.0e-03*exp((((-y[0])-9.0e+01)/8.0e+01))));
}

double Cell::alpha_n (double t, double *y)
{
	return (((1.0e-04*((-y[0])-5.0e+01))/(exp((((-y[0])-5.0e+01)/1.0e+01))-1.0e+00)));
}

double Cell::dndt (int k, double t, double *y)
{
	return ((alpha_n(t,y)*(1.0e+00-y[3]))-(beta_n(t,y)*y[3]));
}
/* ====================================================================================================== */



void setBeats (int nbeats)
{
    __nbeats = nbeats;
}