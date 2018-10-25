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

#include <cstdio>
#include "../include/timer.h"
#include "../include/model.h"

using namespace std;

int main (int argc, char *argv[])
{
  
  if (argc-1 != 1)
  {
    Usage(argv[0]);
    return 1;
  }
  else
  {
    double start, finish, elapsed;

    GET_TIME(start);
    solveModel(argc,argv);
    GET_TIME(finish);
    elapsed = finish - start;
  
    printf("==========================================================\n");
    printf("[!] Time elapsed = %.10lf seconds\n",elapsed);
    printf("==========================================================\n");

    return 0;
  }
}
