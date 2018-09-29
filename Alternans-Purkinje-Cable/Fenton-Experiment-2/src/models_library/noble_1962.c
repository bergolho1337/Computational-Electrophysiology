#include "noble_1962.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    // Original: -87,0.01,0.8,0.01
    // 500ms pacing steady state: -76.1302f, 0.0586201f, 0.741082f, 0.475923f 
    // 450ms pacing steady-state: -78.8607f, 0.050476f, 0.806294f, 0.518928f
    sv[0] = -87.0f;    // V millivolt 
    sv[1] = 0.01f;    // m dimensionless
    sv[2] = 0.8f;    // h dimensionless
    sv[3] = 0.01f;    // n dimensionless
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) {

    uint32_t sv_id;

	int i;

    #pragma omp parallel for private(sv_id)
    for (i = 0; i < num_cells_to_solve; i++) 
    {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i]);

        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current);

    // Solving model using Explicit Euler
    // TO DO: Change to use Rush-Larsen ...
    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current) {

    //State variables
    const real V_old_ = sv[0];
    const real m_old_ = sv[1];
    const real h_old_ = sv[2];
    const real n_old_ = sv[3]; 

    //___________________________________________________________________________
    //Parameters (miliseconds)
    const real Cm = 12.0f;                                 // (microF)
    const real g_na_max = 400.0f;                       // (microS)
    const real E_na = 40.0f;                               // (millivolt)
    const real g_L = 0.075f;                                // (microS)
    const real E_L = -60.0f;                               // (millivolt)

    real calc_I_stim = stim_current;

    // Algebraics
    real g_na =  pow(m_old_, 3.00000)*h_old_*g_na_max;
    real alpha_m = (((1.0e-01*((-V_old_)-4.8e+01))/(exp((((-V_old_)-4.8e+01)/1.5e+01))-1.0e+00)));
    real alpha_h =  ((1.7e-01*exp((((-V_old_)-9.0e+01)/2.0e+01))));
    real alpha_n = (((1.0e-04*((-V_old_)-5.0e+01))/(exp((((-V_old_)-5.0e+01)/1.0e+01))-1.0e+00)));
    real i_na =  (g_na+1.4e-01)*(V_old_ - E_na);
    //real i_na_no_oscilation = (g_na+122.500)*(V_old_ - E_na);
    real beta_m = (((1.2e-01*(V_old_+8.0e+00))/(exp(((V_old_+8.0e+00)/5.0e+00))-1.0e+00)));
    real beta_h = ((1.0/(1.0e+00+exp((((-V_old_)-4.2e+01)/1.0e+01)))));
    real beta_n =  ((2.0e-03*exp((((-V_old_)-9.0e+01)/8.0e+01))));
    //real g_K1 =  1.3f*exp((- V_old_ - 90.0000)/50.0000)+ 0.015f*exp((V_old_+90.0000)/60.0000);
    real g_K1 = (((1.2*exp((((-V_old_)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((V_old_+9.0e+01)/6.0e+01)))));
    real g_K2 =  1.2f*pow(n_old_, 4.00000);
    real i_k =  (g_K1+g_K2)*(V_old_+100.000);
    real i_leak =  g_L*(V_old_ - E_L);

    // Rates
    rDY_[0] = ( - (i_na + i_k + i_leak + calc_I_stim)) / Cm;
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  (beta_m*m_old_) );
    rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  (beta_h*h_old_) );
    rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  (beta_n*n_old_) );

}