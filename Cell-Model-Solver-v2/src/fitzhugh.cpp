#include "../include/fitzhugh.h"

void FitzHugh::initConst (double CONSTANTS[], double RATES[], double STATES[])
{
    STATES[0] = 0;
    STATES[1] = 0;
    // Steady state initial condition
    //STATES[0] = -0.38772821;
    //STATES[1] = 0.15480612;
    CONSTANTS[0] = -0.1;
    CONSTANTS[1] = 3;
    CONSTANTS[2] = 0.005;
}

void FitzHugh::compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{
    RATES[1] =  1.00000*CONSTANTS[2]*(STATES[0] -  CONSTANTS[1]*STATES[1]);
    RATES[0] =  1.00000*(( STATES[0]*(STATES[0] - CONSTANTS[0])*(1.00000 - STATES[0]) - STATES[1])+ALGEBRAIC[0]);
}
void FitzHugh::compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{
    // External stimulus equals: I = 1.0000
    ALGEBRAIC[0] = (t>=0.00000 && t<=0.500000 ?  1.0000 : 0.00000);
    // No external stimulus
    //ALGEBRAIC[0] = (t>=0.00000 && t<=0.500000 ?  0.0000 : 0.00000);
}

void FitzHugh::compRates_FE (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{

}

void FitzHugh::RL (double DT, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{

}

/*
   There are a total of 1 entries in the algebraic variable array.
   There are a total of 2 entries in each of the rate and state variable arrays.
   There are a total of 3 entries in the constant variable array.
*/

/*
 * VOI is t in component Main (millisecond).
 * STATES[0] is v in component Main (dimensionless).
 * STATES[1] is w in component Main (dimensionless).
 * CONSTANTS[0] is alpha in component Main (dimensionless).
 * CONSTANTS[1] is gamma in component Main (dimensionless).
 * CONSTANTS[2] is epsilon in component Main (dimensionless).
 * ALGEBRAIC[0] is I in component Main (dimensionless).
 * RATES[0] is d/dt v in component Main (dimensionless).
 * RATES[1] is d/dt w in component Main (dimensionless).
 */