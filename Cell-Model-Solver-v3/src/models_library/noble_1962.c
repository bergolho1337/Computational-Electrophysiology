#include "noble_1962.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    sv[0] = -75.5344986658f;    // V millivolt 
    sv[1] = 0.060546727200f;    // m dimensionless
    sv[2] = 0.725900135500f;    // h dimensionless
    sv[3] = 0.470923970800f;    // n dimensionless
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) {

    for (int j = 0; j < num_steps; ++j) 
    {
        solve_model_ode_cpu(dt, sv, stim_currents);
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current);

    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current) {

    //State variables
    const real V_old_ = sv[0];
    const real m_old_ = sv[1];
    const real h_old_ = sv[2];
    const real n_old_ = sv[3]; 

    /*
    //Parameters (ms)
    const real Cm = 12.00000000000000000e+00f;             // (microF)
    const real g_na_max = 400000.00000000000000000e+00f;   // (microS)
    const real E_na = 40.00000000000000000e+00f;           // (millivolt)
    const real g_L = 75.00000000000000000e+00f;            // (microS)
    const real E_L = -60.00000000000000000e+00f;           // (millivolt)

    real calc_I_stim = stim_current;

    // Algebraics
    real g_na =  pow(m_old_, 3.00000)*h_old_*g_na_max;
    real alpha_m = ( 100.000*(- V_old_ - 48.0000))/(exp((- V_old_ - 48.0000)/15.0000) - 1.00000);
    real alpha_h =  170.000*exp((- V_old_ - 90.0000)/20.0000);
    real alpha_n = ( 0.100000*(- V_old_ - 50.0000))/(exp((- V_old_ - 50.0000)/10.0000) - 1.00000);
    real i_na =  (g_na+140.000)*(V_old_ - E_na);
    //real i_na_no_oscilation = (g_na+122.500)*(V_old_ - E_na);
    real beta_m = ( 120.000*(V_old_+8.00000))/(exp((V_old_+8.00000)/5.00000) - 1.00000);
    real beta_h = 1000.00/(1.00000+exp((- V_old_ - 42.0000)/10.0000));
    real beta_n =  2.00000*exp((- V_old_ - 90.0000)/80.0000);
    real g_K1 =  1300.00*exp((- V_old_ - 90.0000)/50.0000)+ 15.0000*exp((V_old_+90.0000)/60.0000);
    real g_K2 =  1200.00*pow(n_old_, 4.00000);
    real i_k =  (g_K1+g_K2)*(V_old_+100.000);
    real i_leak =  g_L*(V_old_ - E_L);

    // Rates
    rDY_[0] = (- (i_na + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  beta_m*m_old_) * 1.0E-03;
    rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  beta_h*h_old_) * 1.0E-03;
    rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  beta_n*n_old_) * 1.0E-03;
    */

    //___________________________________________________________________________
    //Parameters (seconds)
    const real Cm = 12.0f;                                 // (microF)
    const real g_na_max = 400000.0f;                       // (microS)
    const real E_na = 40.0f;                               // (millivolt)
    const real g_L = 75.0f;                                // (microS)
    const real E_L = -60.0f;                               // (millivolt)

    real calc_I_stim = stim_current;

    // Algebraics
    real g_na =  pow(m_old_, 3.00000)*h_old_*g_na_max;
    real alpha_m = ( 100.0f*(- V_old_ - 48.0000))/(exp((- V_old_ - 48.0000)/15.0000) - 1.00000);
    real alpha_h =  170.0f*exp((- V_old_ - 90.0000)/20.0000);
    real alpha_n = ( 0.0001*(- V_old_ - 50.0000))/(exp((- V_old_ - 50.0000)/10.0000) - 1.00000);
    real i_na =  (g_na+140.0f)*(V_old_ - E_na);
    //real i_na_no_oscilation = (g_na+122.500)*(V_old_ - E_na);
    real beta_m = ( 120.0f*(V_old_+8.00000))/(exp((V_old_+8.00000)/5.00000) - 1.00000);
    real beta_h = 1000.0f/(1.00000+exp((- V_old_ - 42.0000)/10.0000));
    real beta_n =  0.002*exp((- V_old_ - 90.0000)/80.0000);
    real g_K1 =  1200.0f*exp((- V_old_ - 90.0000)/50.0000)+ 15.0f*exp((V_old_+90.0000)/60.0000);
    real g_K2 =  1200.0f*pow(n_old_, 4.00000);
    real i_k =  (g_K1+g_K2)*(V_old_+100.000);
    real i_leak =  g_L*(V_old_ - E_L);

    // Rates
    rDY_[0] = (- (i_na + i_k + i_leak + calc_I_stim)/Cm);
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  beta_m*m_old_);
    rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  beta_h*h_old_);
    rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  beta_n*n_old_);

}
