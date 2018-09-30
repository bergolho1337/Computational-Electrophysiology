/*
========= FINITE VOLUME METHOD ----- MONODOMAIN EQUATION (CABLE) =============================================
  Problem: solve the cable equation using the monodomain equation applied to the FVM
    { BETA*Cm*V_t = SIGMA*V_xx
    { V'(0,t) = V'(1,t) = I_PMJ           Non-homogeneus Neumann condition at the terminal volumes
    { V(x,0) = V_inf
***********************************************************************************************************
    BETA = Ratio surface-area per volume of the cell (cm^-1)
    Cm = Membrane capacitance of the celullar membrane (uF/cm^2)
    SIGMA = Conductivity of the celullar membrane (mS/cm^2) -- All the fiber has the same conductivity
***********************************************************************************************************
==============================================================================================================
*/

#include <iostream>
#include "monodomain/monodomain_solver.h"
#include "monodomain/ode_solver.h"

using namespace std;

int main (int argc, char *argv[])
{
  struct user_options *options;
  options = new_user_options ();

  struct grid *grid;
  grid = new_grid ();

  struct monodomain_solver *monodomain_solver;
  monodomain_solver = new_monodomain_solver ();

  struct ode_solver *ode_solver;
  ode_solver = new_ode_solver ();

  cout << "Hello world" << endl;
  return 0;
}
