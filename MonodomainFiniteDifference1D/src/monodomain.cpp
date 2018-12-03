#include "../include/monodomain.h"

void solve_diffusion (double *sv, const double dx, const double dt,\
 			const int np, const int nodes,\
			double *vm)
{
	static const double D = 2.5e-04;
	static const double Cm = 12.0;
	static const double r = (dt*D)/(dx*dx);
	
	#pragma omp parallel 
	for (int i = 0; i < np; i++)
	{
		// Boundary node
		if (i == 0)
			vm[i] = (1.0 - r)*sv[i*nodes] + (r)*sv[(i+1)*nodes];
		// Boundary node		
		else if (i == np-1)
			vm[i] = (1.0 - r)*sv[i*nodes] + (r)*sv[(i-1)*nodes];
		// Interior node
		else
			vm[i] = (1.0 - 2.0*r)*sv[i*nodes] + (r)*sv[(i+1)*nodes] + (r)*sv[(i-1)*nodes];
	} 
	
}

void update_state_vector (double *sv, const double *vm,\
			  const int np, const int nodes)
{
	#pragma omp parallel
	for (int i = 0; i < np; i++)
	{
		sv[i*nodes] = vm[i];
	}
}

void solve_reaction (double *sv, double *stims, const double t,\
		     const int np, const int nodes, const double dt)
{
	#pragma omp parallel
	for (int i = 0; i < np; i++)
	{
		double V_old = sv[i*nodes];
		double m_old = sv[i*nodes+1];
		double h_old = sv[i*nodes+2];
		double n_old = sv[i*nodes+3];

		sv[i*nodes] = V_old + dt*dvdt(V_old,m_old,h_old,n_old,stims[i]);
		sv[i*nodes+1] = m_old + dt*dmdt(V_old,m_old);
		sv[i*nodes+2] = h_old + dt*dhdt(V_old,h_old);
		sv[i*nodes+3] = n_old + dt*dndt(V_old,n_old);
		
	}
}
