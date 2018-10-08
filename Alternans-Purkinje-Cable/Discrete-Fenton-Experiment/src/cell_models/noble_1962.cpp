#include "noble_1962.h"

void init_cell_model_data (uint32_t *num_equations) 
{
    *num_equations = NEQ;
}

void set_model_initial_conditions_cpu (double *y_old)
{
    // Original: -75.5344986658,0.0605467272,0.7259001355,0.4709239708
    y_old[0] = -75.5344986658;    // V millivolt 
    y_old[1] = 0.0605467272;      // m dimensionless
    y_old[2] = 0.7259001355;      // h dimensionless
    y_old[3] = 0.4709239708;      // n dimensionless
}

void solve_model_odes_cpu (const double dt, double *sv, double *vstar,\
                            double *stim_currents, const uint32_t num_cells_to_solve) 
{
    uint32_t sv_id;

	uint32_t i;

    for (i = 0; i < num_cells_to_solve; i++) 
    {
        sv_id = i;
        solve_model_ode_cpu(dt, sv + (sv_id * NEQ), vstar[i], stim_currents[i]);
    }
}

void solve_model_ode_cpu(double dt, double *sv, double vstar, double stim_current)
{
    double rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, vstar, rDY, stim_current);

    // Solving model using Explicit Euler
    sv[0] = dt*rDY[0] + vstar;
    for(int i = 1; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void RHS_cpu(const double *sv, const double vstar, double *rDY_, double stim_current)
{
    //State variables
    const double V_old_ = sv[0];
    const double m_old_ = sv[1];
    const double h_old_ = sv[2];
    const double n_old_ = sv[3];
    const double V_star_ = vstar; 

    //___________________________________________________________________________
    //Parameters (miliseconds)
    const double Cm = 12.0;                                 // (microF)
    const double g_na_max = 400.0;                       // (microS)
    const double E_na = 40.0;                               // (millivolt)
    const double g_L = 0.075;                                // (microS)
    const double E_L = -60.0;                               // (millivolt)

    double calc_I_stim = stim_current;

    // Algebraics
    double g_na =  pow(m_old_, 3.00000)*h_old_*g_na_max;
    double alpha_m = (((1.0e-01*((-V_old_)-4.8e+01))/(exp((((-V_old_)-4.8e+01)/1.5e+01))-1.0e+00)));
    double alpha_h =  ((1.7e-01*exp((((-V_old_)-9.0e+01)/2.0e+01))));
    double alpha_n = (((1.0e-04*((-V_old_)-5.0e+01))/(exp((((-V_old_)-5.0e+01)/1.0e+01))-1.0e+00)));
    double i_na =  (g_na+1.4e-01)*(V_star_ - E_na);
    //double i_na_no_oscilation = (g_na+122.500)*(V_old_ - E_na);
    double beta_m = (((1.2e-01*(V_old_+8.0e+00))/(exp(((V_old_+8.0e+00)/5.0e+00))-1.0e+00)));
    double beta_h = ((1.0/(1.0e+00+exp((((-V_old_)-4.2e+01)/1.0e+01)))));
    double beta_n =  ((2.0e-03*exp((((-V_old_)-9.0e+01)/8.0e+01))));
    //double g_K1 =  1.3f*exp((- V_old_ - 90.0000)/50.0000)+ 0.015f*exp((V_old_+90.0000)/60.0000);
    double g_K1 = (((1.2*exp((((-V_star_)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((V_star_+9.0e+01)/6.0e+01)))));
    double g_K2 =  1.2*pow(n_old_, 4.00000);
    double i_k =  (g_K1+g_K2)*(V_star_+100.000);
    double i_leak =  g_L*(V_star_ - E_L);

    // Rates
    rDY_[0] = ( (-(i_na + i_k + i_leak) + calc_I_stim ) ) / Cm;
    rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  (beta_m*m_old_) );
    rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  (beta_h*h_old_) );
    rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  (beta_n*n_old_) );
}