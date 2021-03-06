#include "../include/noble.h"

void Noble::initConst (double CONSTANTS[], double RATES[], double STATES[])
{
    STATES[0] = -87;
    STATES[1] = 0.01;
    STATES[2] = 0.8;
    STATES[3] = 0.01;
    CONSTANTS[0] = 12;
    CONSTANTS[1] = 400000;
    CONSTANTS[2] = 40;
    CONSTANTS[3] = 75;
    CONSTANTS[4] = -60;
}

void Noble::compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{
    RATES[0] = (- (ALGEBRAIC[4]+ALGEBRAIC[10]+ALGEBRAIC[11])/CONSTANTS[0]) * 1.0E-03;
    RATES[1] =  (ALGEBRAIC[1]*(1.00000 - STATES[1]) -  ALGEBRAIC[5]*STATES[1]) * 1.0E-03;
    RATES[2] =  (ALGEBRAIC[2]*(1.00000 - STATES[2]) -  ALGEBRAIC[6]*STATES[2]) * 1.0E-03;
    RATES[3] =  (ALGEBRAIC[3]*(1.00000 - STATES[3]) -  ALGEBRAIC[7]*STATES[3]) * 1.0E-03;
}

void Noble::compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{
    ALGEBRAIC[0] =  pow(STATES[1], 3.00000)*STATES[2]*CONSTANTS[1];
    ALGEBRAIC[1] = ( 100.000*(- STATES[0] - 48.0000))/(exp((- STATES[0] - 48.0000)/15.0000) - 1.00000);
    ALGEBRAIC[2] =  170.000*exp((- STATES[0] - 90.0000)/20.0000);
    ALGEBRAIC[3] = ( 0.100000*(- STATES[0] - 50.0000))/(exp((- STATES[0] - 50.0000)/10.0000) - 1.00000);
    ALGEBRAIC[4] =  (ALGEBRAIC[0]+140.000)*(STATES[0] - CONSTANTS[2]);
    ALGEBRAIC[5] = ( 120.000*(STATES[0]+8.00000))/(exp((STATES[0]+8.00000)/5.00000) - 1.00000);
    ALGEBRAIC[6] = 1000.00/(1.00000+exp((- STATES[0] - 42.0000)/10.0000));
    ALGEBRAIC[7] =  2.00000*exp((- STATES[0] - 90.0000)/80.0000);
    ALGEBRAIC[8] =  1200.00*exp((- STATES[0] - 90.0000)/50.0000)+ 15.0000*exp((STATES[0]+90.0000)/60.0000);
    ALGEBRAIC[9] =  1200.00*pow(STATES[3], 4.00000);
    ALGEBRAIC[10] =  (ALGEBRAIC[8]+ALGEBRAIC[9])*(STATES[0]+100.000);
    ALGEBRAIC[11] =  CONSTANTS[3]*(STATES[0] - CONSTANTS[4]);
}

void Noble::compRates_FE (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{
    // V
    RATES[0] = (- (ALGEBRAIC[4]+ALGEBRAIC[10]+ALGEBRAIC[11])/CONSTANTS[0]) * 1.0E-03;
}

void Noble::RL (double DT, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{
    // Forward Euler for V
    STATES[0] = STATES[0] + RATES[0]*DT;
    
    // Rush-Larsen for m, h, n
    // m_inf = alpha_m / (alpha_m + beta_m)
    double m_inf = ALGEBRAIC[1] / (ALGEBRAIC[1] + ALGEBRAIC[5]);
    // tau_m = 1.0 / (alpha_m + beta_m)
    double tau_m = 1.0 / (ALGEBRAIC[1] + ALGEBRAIC[5]);
    STATES[1] = m_inf + ((STATES[1] - m_inf)*exp(-DT * 1.0E-03/tau_m));

    // h_inf = alpha_h / (alpha_h + beta_h)
    double h_inf = ALGEBRAIC[2] / (ALGEBRAIC[2] + ALGEBRAIC[6]);
    // tau_h = 1.0 / (alpha_h + beta_h)
    double tau_h = 1.0 / (ALGEBRAIC[2] + ALGEBRAIC[6]);
    STATES[2] = h_inf + ((STATES[2] - h_inf)*exp(-DT * 1.0E-03/tau_h));

    // n_inf = alpha_n / (alpha_n + beta_n)
    double n_inf = ALGEBRAIC[3] / (ALGEBRAIC[3] + ALGEBRAIC[7]);
    // tau_n = 1.0 / (alpha_n + beta_n)
    double tau_n = 1.0 / (ALGEBRAIC[3] + ALGEBRAIC[7]);
    STATES[3] = n_inf + ((STATES[3] - n_inf)*exp(-DT * 1.0E-03/tau_n));
}

/*
   There are a total of 12 entries in the algebraic variable array.
   There are a total of 4 entries in each of the rate and state variable arrays.
   There are a total of 5 entries in the constant variable array.
*/

/*
 * VOI is time in component environment (milisecond).
 * STATES[0] is V in component membrane (millivolt).
 * ALGEBRAIC[0] is g_Na in component sodium_channel (microS).
 * ALGEBRAIC[1] is alpha_m in component sodium_channel_m_gate (per_milisecond).
 * ALGEBRAIC[2] is alpha_h in component sodium_channel_h_gate (per_milisecond).
 * ALGEBRAIC[3] is alpha_n in component potassium_channel_n_gate (per_milisecond).
 * ALGEBRAIC[4] is i_Na in component sodium_channel (nanoA).
 * ALGEBRAIC[5] is beta_m in component sodium_channel_m_gate (per_milisecond).
 * ALGEBRAIC[6] is beta_h in component sodium_channel_h_gate (per_milisecond).
 * ALGEBRAIC[7] is beta_n in component potassium_channel_n_gate (per_milisecond).
 * ALGEBRAIC[8] is g_K1 in component potassium_channel (microS).
 * ALGEBRAIC[9] is g_K2 in component potassium_channel (microS).
 * ALGEBRAIC[10] is i_K in component potassium_channel (nanoA).
 * ALGEBRAIC[11] is i_Leak in component leakage_current (nanoA).
 * STATES[1] is m in component sodium_channel_m_gate (dimensionless).
 * STATES[2] is h in component sodium_channel_h_gate (dimensionless).
 * STATES[3] is n in component potassium_channel_n_gate (dimensionless).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[1] is d/dt m in component sodium_channel_m_gate (dimensionless).
 * RATES[2] is d/dt h in component sodium_channel_h_gate (dimensionless).
 * RATES[3] is d/dt n in component potassium_channel_n_gate (dimensionless).
 * CONSTANTS[0] is Cm in component membrane (microF).
 * CONSTANTS[1] is g_Na_max in component sodium_channel (microS).
 * CONSTANTS[2] is E_Na in component sodium_channel (millivolt).
 * CONSTANTS[3] is g_L in component leakage_current (microS).
 * CONSTANTS[4] is E_L in component leakage_current (millivolt).
 * 
 */
