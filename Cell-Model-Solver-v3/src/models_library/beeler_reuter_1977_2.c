/*
   There are a total of 18 entries in the algebraic variable array.
   There are a total of 8 entries in each of the rate and state variable arrays.
   There are a total of 10 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (ms).
 * STATES[0] is V in component membrane (mV).
 * STATES[1] is m in component sodium_current_m_gate (dimensionless).
 * STATES[2] is h in component sodium_current_h_gate (dimensionless).
 * STATES[3] is j in component sodium_current_j_gate (dimensionless).
 * STATES[4] is Cai in component slow_inward_current (concentration_units).
 * STATES[5] is d in component slow_inward_current_d_gate (dimensionless).
 * STATES[6] is f in component slow_inward_current_f_gate (dimensionless).
 * STATES[7] is x1 in component time_dependent_outward_current_x1_gate (dimensionless).
 
 * CONSTANTS[0] is C in component membrane (uF_per_mm2).
 * CONSTANTS[1] is g_Na in component sodium_current (mS_per_mm2).
 * CONSTANTS[2] is E_Na in component sodium_current (mV).
 * CONSTANTS[3] is g_Nac in component sodium_current (mS_per_mm2).
 * CONSTANTS[4] is g_s in component slow_inward_current (mS_per_mm2).
 * CONSTANTS[5] is IstimStart in component stimulus_protocol (ms).
 * CONSTANTS[6] is IstimEnd in component stimulus_protocol (ms).
 * CONSTANTS[7] is IstimAmplitude in component stimulus_protocol (uA_per_mm2).
 * CONSTANTS[8] is IstimPeriod in component stimulus_protocol (ms).
 * CONSTANTS[9] is IstimPulseDuration in component stimulus_protocol (ms).
 
 * ALGEBRAIC[0] is i_Na in component sodium_current (uA_per_mm2).
 * ALGEBRAIC[1] is alpha_m in component sodium_current_m_gate (per_ms).
 * ALGEBRAIC[2] is alpha_h in component sodium_current_h_gate (per_ms).
 * ALGEBRAIC[3] is alpha_j in component sodium_current_j_gate (per_ms).
 * ALGEBRAIC[4] is alpha_d in component slow_inward_current_d_gate (per_ms).
 * ALGEBRAIC[5] is alpha_f in component slow_inward_current_f_gate (per_ms).
 * ALGEBRAIC[6] is alpha_x1 in component time_dependent_outward_current_x1_gate (per_ms).
 * ALGEBRAIC[7] is E_s in component slow_inward_current (mV).
 * ALGEBRAIC[8] is beta_m in component sodium_current_m_gate (per_ms).
 * ALGEBRAIC[9] is beta_h in component sodium_current_h_gate (per_ms).
 * ALGEBRAIC[10] is beta_j in component sodium_current_j_gate (per_ms).
 * ALGEBRAIC[11] is beta_d in component slow_inward_current_d_gate (per_ms).
 * ALGEBRAIC[12] is beta_f in component slow_inward_current_f_gate (per_ms).
 * ALGEBRAIC[13] is beta_x1 in component time_dependent_outward_current_x1_gate (per_ms).
 * ALGEBRAIC[14] is i_s in component slow_inward_current (uA_per_mm2).
 * ALGEBRAIC[15] is i_x1 in component time_dependent_outward_current (uA_per_mm2).
 * ALGEBRAIC[16] is i_K1 in component time_independent_outward_current (uA_per_mm2).
 * ALGEBRAIC[17] is Istim in component stimulus_protocol (uA_per_mm2).
 
 * RATES[0] is d/dt V in component membrane (mV).
 * RATES[1] is d/dt m in component sodium_current_m_gate (dimensionless).
 * RATES[2] is d/dt h in component sodium_current_h_gate (dimensionless).
 * RATES[3] is d/dt j in component sodium_current_j_gate (dimensionless).
 * RATES[4] is d/dt Cai in component slow_inward_current (concentration_units).
 * RATES[5] is d/dt d in component slow_inward_current_d_gate (dimensionless).
 * RATES[6] is d/dt f in component slow_inward_current_f_gate (dimensionless).
 * RATES[7] is d/dt x1 in component time_dependent_outward_current_x1_gate (dimensionless).
 */
void
initConsts(double* CONSTANTS, double* RATES, double *STATES)
{
STATES[0] = -84.624;
STATES[1] = 0.011;
STATES[2] = 0.988;
STATES[3] = 0.975;
STATES[4] = 1e-4;
STATES[5] = 0.003;
STATES[6] = 0.994;
STATES[7] = 0.0001;

CONSTANTS[0] = 0.01;
CONSTANTS[1] = 4e-2;
CONSTANTS[2] = 50;
CONSTANTS[3] = 3e-5;
CONSTANTS[4] = 9e-4;
CONSTANTS[5] = 10;
CONSTANTS[6] = 50000;
CONSTANTS[7] = 0.5;
CONSTANTS[8] = 1000;
CONSTANTS[9] = 1;
}
void
computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
RATES[1] =  ALGEBRAIC[1]*(1.00000 - STATES[1]) -  ALGEBRAIC[8]*STATES[1];
RATES[2] =  ALGEBRAIC[2]*(1.00000 - STATES[2]) -  ALGEBRAIC[9]*STATES[2];
RATES[3] =  ALGEBRAIC[3]*(1.00000 - STATES[3]) -  ALGEBRAIC[10]*STATES[3];
RATES[5] =  ALGEBRAIC[4]*(1.00000 - STATES[5]) -  ALGEBRAIC[11]*STATES[5];
RATES[6] =  ALGEBRAIC[5]*(1.00000 - STATES[6]) -  ALGEBRAIC[12]*STATES[6];
RATES[7] =  ALGEBRAIC[6]*(1.00000 - STATES[7]) -  ALGEBRAIC[13]*STATES[7];
RATES[4] = ( - 0.0100000*ALGEBRAIC[14])/1.00000+ 0.0700000*(0.000100000 - STATES[4]);
RATES[0] = (ALGEBRAIC[17] - (ALGEBRAIC[0]+ALGEBRAIC[14]+ALGEBRAIC[15]+ALGEBRAIC[16]))/CONSTANTS[0];
}
void
computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
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

#include "beeler_reuter_1977.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    sv[0] = -84.624;        // V (mV)
    sv[1] = 0.011;          // m (dimensionless)
    sv[2] = 0.988;          // h (dimensionless)
    sv[3] = 0.975;          // j (dimensionless)
    sv[4] = 1e-4;           // Cai (concentration_units)
    sv[5] = 0.003;          // d (dimensionless)
    sv[6] = 0.994;          // f (dimensionless)
    sv[7] = 0.0001;         // x1 (dimensionless)
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) {

    for (int j = 0; j < num_steps; ++j) 
    {
        solve_model_ode_cpu(dt, sv, stim_currents);
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    real rY[NEQ], rDY[NEQ];

    // Save old value of the state vector
    for(int i = 0; i < NEQ; i++)
    {
        rY[i] = sv[i];
    }
        

    // Compute Right-hand-side of the ODE's
    RHS_cpu(rY, rDY, stim_current);

    // Solve model using Forward Euler
    for(int i = 0; i < NEQ; i++)
    {
        sv[i] = dt*rDY[i] + rY[i];
    }
        
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current) {

    // Constants
    const real C = 0.01;        // Component membrane (uF_per_mm2)
    const real g_na = 4e-2;     // Component sodium current (mS_per_mm2)
    const real E_na = 50;       // Component sodium current (mV)
    const real g_nac = 3e-5;    // Component sodium current (mS_per_mm2)
    const real g_s = 9e-4;      // Component slow inward current (mS_per_mm2)

    //State variables
    const real V_old_ = -84.624;        // V (mV)
    const real m_old_ = 0.011;          // m (dimensionless)
    const real h_old_ = 0.988;          // h (dimensionless)
    const real j_old_ = 0.975;          // j (dimensionless)
    const real Cai_old_ = 1e-4;           // Cai (concentration_units)
    const real d_old_ = 0.003;          // d (dimensionless)
    const real f_old_ = 0.994;          // f (dimensionless)
    const real x1_old_ = 0.0001;         // x1 (dimensionless)

    // Algebraics
    real alpha_m = ( - 1.00000*(V_old_+47.0000))/(exp( - 0.100000*(V_old_+47.0000)) - 1.00000);
    real beta_m =  40.0000*exp( - 0.0560000*(V_old_+72.0000));
    real alpha_h =  0.126000*exp( - 0.250000*(V_old_+77.0000));
    real beta_h = 1.70000/(exp( - 0.0820000*(V_old_+22.5000))+1.00000);
    real alpha_j = ( 0.0550000*exp( - 0.250000*(V_old_+78.0000)))/(exp( - 0.200000*(V_old_+78.0000))+1.00000);
    real beta_j = 0.300000/(exp( - 0.100000*(V_old_+32.0000))+1.00000);
    real alpha_d = ( 0.0950000*exp(- (V_old_ - 5.00000)/100.000))/(1.00000+exp(- (V_old_ - 5.00000)/13.8900));
    real beta_d = ( 0.0700000*exp(- (V_old_+44.0000)/59.0000))/(1.00000+exp((V_old_+44.0000)/20.0000));
    real alpha_f = ( 0.0120000*exp(- (V_old_+28.0000)/125.000))/(1.00000+exp((V_old_+28.0000)/6.67000));
    real beta_f = ( 0.00650000*exp(- (V_old_+30.0000)/50.0000))/(1.00000+exp(- (V_old_+30.0000)/5.00000));
    real alpha_x1 = ( 0.000500000*exp((V_old_+50.0000)/12.1000))/(1.00000+exp((V_old_+50.0000)/17.5000));
    real beta_x1 = ( 0.00130000*exp(- (V_old_+20.0000)/16.6700))/(1.00000+exp(- (V_old_+20.0000)/25.0000));
    real E_s = - 82.3000 -  13.0287*log( Cai_old_*0.00100000);
    real i_s =  g_s*d_old_*f_old_*(V_old_ - E_s);
    real i_na =  ( g_na*pow(m_old_, 3.00000)*h_old_*j_old_+g_nac)*(V_old_ - E_na);
    real i_x1 = ( x1_old_*0.00800000*(exp( 0.0400000*(V_old_+77.0000)) - 1.00000))/exp( 0.0400000*(V_old_+35.0000));
    real i_k1 =  0.00350000*(( 4.00000*(exp( 0.0400000*(V_old_+85.0000)) - 1.00000))/(exp( 0.0800000*(V_old_+53.0000))+exp( 0.0400000*(V_old_+53.0000)))+( 0.200000*(V_old_+23.0000))/(1.00000 - exp( - 0.0400000*(V_old_+23.0000))));
    real i_stim = stim_current;
    
    // Rates
    real d_dt_V = (i_stim - (i_na+i_s+i_x1+i_k1))/C;
    real d_dt_m =  alpha_m*(1.00000 - m_old_) -  beta_m*m_old_;
    real d_dt_h =  alpha_h*(1.00000 - h_old_) -  beta_h*h_old_;
    RATES[3] =  ALGEBRAIC[3]*(1.00000 - STATES[3]) -  ALGEBRAIC[10]*STATES[3];
    RATES[5] =  ALGEBRAIC[4]*(1.00000 - STATES[5]) -  ALGEBRAIC[11]*STATES[5];
    RATES[6] =  ALGEBRAIC[5]*(1.00000 - STATES[6]) -  ALGEBRAIC[12]*STATES[6];
    RATES[7] =  ALGEBRAIC[6]*(1.00000 - STATES[7]) -  ALGEBRAIC[13]*STATES[7];
    RATES[4] = ( - 0.0100000*ALGEBRAIC[14])/1.00000+ 0.0700000*(0.000100000 - STATES[4]);
    rDY_[0] = Jstim + Jin + Jout;
    rDY_[1] = (V_old_<V_gate ? (1.00000 - h_old_)/tau_open : - h_old_/tau_close);

}