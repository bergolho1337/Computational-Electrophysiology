/*
 * VOI is time in component environment (ms).
 * STATES[0] is V in component membrane (mV).
 * CONSTANTS[0] is C in component membrane (uF_per_mm2).
 * ALGEBRAIC[0] is i_Na in component sodium_current (uA_per_mm2).
 * ALGEBRAIC[14] is i_s in component slow_inward_current (uA_per_mm2).
 * ALGEBRAIC[15] is i_x1 in component time_dependent_outward_current (uA_per_mm2).
 * ALGEBRAIC[16] is i_K1 in component time_independent_outward_current (uA_per_mm2).
 * ALGEBRAIC[17] is Istim in component stimulus_protocol (uA_per_mm2).
 * CONSTANTS[1] is g_Na in component sodium_current (mS_per_mm2).
 * CONSTANTS[2] is E_Na in component sodium_current (mV).
 * CONSTANTS[3] is g_Nac in component sodium_current (mS_per_mm2).
 * STATES[1] is m in component sodium_current_m_gate (dimensionless).
 * STATES[2] is h in component sodium_current_h_gate (dimensionless).
 * STATES[3] is j in component sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[1] is alpha_m in component sodium_current_m_gate (per_ms).
 * ALGEBRAIC[8] is beta_m in component sodium_current_m_gate (per_ms).
 * ALGEBRAIC[2] is alpha_h in component sodium_current_h_gate (per_ms).
 * ALGEBRAIC[9] is beta_h in component sodium_current_h_gate (per_ms).
 * ALGEBRAIC[3] is alpha_j in component sodium_current_j_gate (per_ms).
 * ALGEBRAIC[10] is beta_j in component sodium_current_j_gate (per_ms).
 * CONSTANTS[4] is g_s in component slow_inward_current (mS_per_mm2).
 * ALGEBRAIC[7] is E_s in component slow_inward_current (mV).
 * STATES[4] is Cai in component slow_inward_current (concentration_units).
 * STATES[5] is d in component slow_inward_current_d_gate (dimensionless).
 * STATES[6] is f in component slow_inward_current_f_gate (dimensionless).
 * ALGEBRAIC[4] is alpha_d in component slow_inward_current_d_gate (per_ms).
 * ALGEBRAIC[11] is beta_d in component slow_inward_current_d_gate (per_ms).
 * ALGEBRAIC[5] is alpha_f in component slow_inward_current_f_gate (per_ms).
 * ALGEBRAIC[12] is beta_f in component slow_inward_current_f_gate (per_ms).
 * STATES[7] is x1 in component time_dependent_outward_current_x1_gate (dimensionless).
 * ALGEBRAIC[6] is alpha_x1 in component time_dependent_outward_current_x1_gate (per_ms).
 * ALGEBRAIC[13] is beta_x1 in component time_dependent_outward_current_x1_gate (per_ms).
 * CONSTANTS[5] is IstimStart in component stimulus_protocol (ms).
 * CONSTANTS[6] is IstimEnd in component stimulus_protocol (ms).
 * CONSTANTS[7] is IstimAmplitude in component stimulus_protocol (uA_per_mm2).
 * CONSTANTS[8] is IstimPeriod in component stimulus_protocol (ms).
 * CONSTANTS[9] is IstimPulseDuration in component stimulus_protocol (ms).
 * RATES[0] is d/dt V in component membrane (mV).
 * RATES[1] is d/dt m in component sodium_current_m_gate (dimensionless).
 * RATES[2] is d/dt h in component sodium_current_h_gate (dimensionless).
 * RATES[3] is d/dt j in component sodium_current_j_gate (dimensionless).
 * RATES[4] is d/dt Cai in component slow_inward_current (concentration_units).
 * RATES[5] is d/dt d in component slow_inward_current_d_gate (dimensionless).
 * RATES[6] is d/dt f in component slow_inward_current_f_gate (dimensionless).
 * RATES[7] is d/dt x1 in component time_dependent_outward_current_x1_gate (dimensionless).
 */

#include <iostream>
#include <cmath>

using namespace std;

/*
   There are a total of 18 entries in the algebraic variable array.
   There are a total of 8 entries in each of the rate and state variable arrays.
   There are a total of 10 entries in the constant variable array.
 */
const int NRATES = 8;
const int NALGEBRAIC = 18;
const int NCONSTANTS = 10;

void initConsts(double* CONSTANTS, double* RATES, double *STATES)
{
    STATES[0] = -84.624;
    CONSTANTS[0] = 0.01;
    CONSTANTS[1] = 4e-2;
    CONSTANTS[2] = 50;
    CONSTANTS[3] = 3e-5;
    STATES[1] = 0.011;
    STATES[2] = 0.988;
    STATES[3] = 0.975;
    CONSTANTS[4] = 9e-4;
    STATES[4] = 1e-4;
    STATES[5] = 0.003;
    STATES[6] = 0.994;
    STATES[7] = 0.0001;
    CONSTANTS[5] = 10;
    CONSTANTS[6] = 50000;
    CONSTANTS[7] = 0.5;
    CONSTANTS[8] = 1000;
    CONSTANTS[9] = 1;
}
void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
    ALGEBRAIC[1] = ( - 1.00000*(STATES[0]+47.0000))/(exp( - 0.100000*(STATES[0]+47.0000)) - 1.00000);
    ALGEBRAIC[8] =  40.0000*exp( - 0.0560000*(STATES[0]+72.0000));
    RATES[1] =  ALGEBRAIC[1]*(1.00000 - STATES[1]) -  ALGEBRAIC[8]*STATES[1];
    ALGEBRAIC[2] =  0.126000*exp( - 0.250000*(STATES[0]+77.0000));
    ALGEBRAIC[9] = 1.70000/(exp( - 0.0820000*(STATES[0]+22.5000))+1.00000);
    RATES[2] =  ALGEBRAIC[2]*(1.00000 - STATES[2]) -  ALGEBRAIC[9]*STATES[2];
    ALGEBRAIC[3] = ( 0.0550000*exp( - 0.250000*(STATES[0]+78.0000)))/(exp( - 0.200000*(STATES[0]+78.0000))+1.00000);
    ALGEBRAIC[10] = 0.300000/(exp( - 0.100000*(STATES[0]+32.0000))+1.00000);
    RATES[3] =  ALGEBRAIC[3]*(1.00000 - STATES[3]) -  ALGEBRAIC[10]*STATES[3];
    ALGEBRAIC[4] = ( 0.0950000*exp(- (STATES[0] - 5.00000)/100.000))/(1.00000+exp(- (STATES[0] - 5.00000)/13.8900));
    ALGEBRAIC[11] = ( 0.0700000*exp(- (STATES[0]+44.0000)/59.0000))/(1.00000+exp((STATES[0]+44.0000)/20.0000));
    RATES[5] =  ALGEBRAIC[4]*(1.00000 - STATES[5]) -  ALGEBRAIC[11]*STATES[5];
    ALGEBRAIC[5] = ( 0.0120000*exp(- (STATES[0]+28.0000)/125.000))/(1.00000+exp((STATES[0]+28.0000)/6.67000));
    ALGEBRAIC[12] = ( 0.00650000*exp(- (STATES[0]+30.0000)/50.0000))/(1.00000+exp(- (STATES[0]+30.0000)/5.00000));
    RATES[6] =  ALGEBRAIC[5]*(1.00000 - STATES[6]) -  ALGEBRAIC[12]*STATES[6];
    ALGEBRAIC[6] = ( 0.000500000*exp((STATES[0]+50.0000)/12.1000))/(1.00000+exp((STATES[0]+50.0000)/17.5000));
    ALGEBRAIC[13] = ( 0.00130000*exp(- (STATES[0]+20.0000)/16.6700))/(1.00000+exp(- (STATES[0]+20.0000)/25.0000));
    RATES[7] =  ALGEBRAIC[6]*(1.00000 - STATES[7]) -  ALGEBRAIC[13]*STATES[7];
    ALGEBRAIC[7] = - 82.3000 -  13.0287*log( STATES[4]*0.00100000);
    ALGEBRAIC[14] =  CONSTANTS[4]*STATES[5]*STATES[6]*(STATES[0] - ALGEBRAIC[7]);
    RATES[4] = ( - 0.0100000*ALGEBRAIC[14])/1.00000+ 0.0700000*(0.000100000 - STATES[4]);
    ALGEBRAIC[0] =  ( CONSTANTS[1]*pow(STATES[1], 3.00000)*STATES[2]*STATES[3]+CONSTANTS[3])*(STATES[0] - CONSTANTS[2]);
    ALGEBRAIC[15] = ( STATES[7]*0.00800000*(exp( 0.0400000*(STATES[0]+77.0000)) - 1.00000))/exp( 0.0400000*(STATES[0]+35.0000));
    ALGEBRAIC[16] =  0.00350000*(( 4.00000*(exp( 0.0400000*(STATES[0]+85.0000)) - 1.00000))/(exp( 0.0800000*(STATES[0]+53.0000))+exp( 0.0400000*(STATES[0]+53.0000)))+( 0.200000*(STATES[0]+23.0000))/(1.00000 - exp( - 0.0400000*(STATES[0]+23.0000))));
    ALGEBRAIC[17] = (VOI>=CONSTANTS[5]&&VOI<=CONSTANTS[6]&&(VOI - CONSTANTS[5]) -  floor((VOI - CONSTANTS[5])/CONSTANTS[8])*CONSTANTS[8]<=CONSTANTS[9] ? CONSTANTS[7] : 0.00000);
    RATES[0] = (ALGEBRAIC[17] - (ALGEBRAIC[0]+ALGEBRAIC[14]+ALGEBRAIC[15]+ALGEBRAIC[16]))/CONSTANTS[0];
}

void computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
    ALGEBRAIC[1] = ( - 1.00000*(STATES[0]+47.0000))/(exp( - 0.100000*(STATES[0]+47.0000)) - 1.00000);
    ALGEBRAIC[8] =  40.0000*exp( - 0.0560000*(STATES[0]+72.0000));
    ALGEBRAIC[2] =  0.126000*exp( - 0.250000*(STATES[0]+77.0000));
    ALGEBRAIC[9] = 1.70000/(exp( - 0.0820000*(STATES[0]+22.5000))+1.00000);
    ALGEBRAIC[3] = ( 0.0550000*exp( - 0.250000*(STATES[0]+78.0000)))/(exp( - 0.200000*(STATES[0]+78.0000))+1.00000);
    ALGEBRAIC[10] = 0.300000/(exp( - 0.100000*(STATES[0]+32.0000))+1.00000);
    ALGEBRAIC[4] = ( 0.0950000*exp(- (STATES[0] - 5.00000)/100.000))/(1.00000+exp(- (STATES[0] - 5.00000)/13.8900));
    ALGEBRAIC[11] = ( 0.0700000*exp(- (STATES[0]+44.0000)/59.0000))/(1.00000+exp((STATES[0]+44.0000)/20.0000));
    ALGEBRAIC[5] = ( 0.0120000*exp(- (STATES[0]+28.0000)/125.000))/(1.00000+exp((STATES[0]+28.0000)/6.67000));
    ALGEBRAIC[12] = ( 0.00650000*exp(- (STATES[0]+30.0000)/50.0000))/(1.00000+exp(- (STATES[0]+30.0000)/5.00000));
    ALGEBRAIC[6] = ( 0.000500000*exp((STATES[0]+50.0000)/12.1000))/(1.00000+exp((STATES[0]+50.0000)/17.5000));
    ALGEBRAIC[13] = ( 0.00130000*exp(- (STATES[0]+20.0000)/16.6700))/(1.00000+exp(- (STATES[0]+20.0000)/25.0000));
    ALGEBRAIC[7] = - 82.3000 -  13.0287*log( STATES[4]*0.00100000);
    ALGEBRAIC[14] =  CONSTANTS[4]*STATES[5]*STATES[6]*(STATES[0] - ALGEBRAIC[7]);
    ALGEBRAIC[0] =  ( CONSTANTS[1]*pow(STATES[1], 3.00000)*STATES[2]*STATES[3]+CONSTANTS[3])*(STATES[0] - CONSTANTS[2]);
    ALGEBRAIC[15] = ( STATES[7]*0.00800000*(exp( 0.0400000*(STATES[0]+77.0000)) - 1.00000))/exp( 0.0400000*(STATES[0]+35.0000));
    ALGEBRAIC[16] =  0.00350000*(( 4.00000*(exp( 0.0400000*(STATES[0]+85.0000)) - 1.00000))/(exp( 0.0800000*(STATES[0]+53.0000))+exp( 0.0400000*(STATES[0]+53.0000)))+( 0.200000*(STATES[0]+23.0000))/(1.00000 - exp( - 0.0400000*(STATES[0]+23.0000))));
    ALGEBRAIC[17] = (VOI>=CONSTANTS[5]&&VOI<=CONSTANTS[6]&&(VOI - CONSTANTS[5]) -  floor((VOI - CONSTANTS[5])/CONSTANTS[8])*CONSTANTS[8]<=CONSTANTS[9] ? CONSTANTS[7] : 0.00000);
}

int main ()
{
    double tmax = 100.0;
    double dt = 0.01;
    int N = nearbyint(tmax / dt);

    double *CONSTANTS = new double[NCONSTANTS];
    double *RATES = new double[NRATES];
    double *STATES = new double[NRATES];
    double *ALGEBRAIC = new double[NALGEBRAIC];

    initConsts(CONSTANTS,RATES,STATES);
    for (int i = 0; i < N; i++)
    {
        double VOI = i*dt;
        computeVariables(VOI,CONSTANTS,RATES,STATES,ALGEBRAIC);
        computeRates(VOI,CONSTANTS,RATES,STATES,ALGEBRAIC);
        cout << "t = " << VOI << " -- V = " << STATES[0] << endl;
    }
    
    delete [] CONSTANTS;
    delete [] RATES;
    delete [] STATES;
    delete [] ALGEBRAIC;

    return 0;
}