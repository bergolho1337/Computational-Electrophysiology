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

void solve_model_odes_cpu (struct cell_data *volumes, double *stim_currents,\
                           const double dt ,const uint32_t num_cells_to_solve) 
{
    uint32_t n_odes = NEQ;
    for (uint32_t i = 0; i < num_cells_to_solve; i++)
    {
        double *y_old = volumes[i].yOld;
        double *y_star = volumes[i].yStar;
        double *y_new = volumes[i].yNew;

        double f1 = dvdt(y_star,stim_currents[i]);
        y_new[0] = y_star[0] + f1*dt;

        double f2 = dmdt(y_old);
        y_new[1] = y_old[1] + f2*dt;

        double f3 = dhdt(y_old);
        y_new[2] = y_old[2] + f3*dt;

        double f4 = dndt(y_old);
        y_new[3] = y_old[3] + f4*dt;
    }
}
// =================================================================================================
// dV / dt
double g_na (const double *y)
{
    const double m_old_ = y[1];
    const double h_old_ = y[2];

    return pow(m_old_, 3.00000)*h_old_*g_na_max;
}

double i_na (const double *y)
{
    const double V_old_ = y[0];

    return (g_na(y) + 1.4e-01)*(V_old_ - E_na);

}

double g_K2 (const double *y)
{
    const double n_old_ = y[3];

    return 1.2*pow(n_old_, 4.00000);
}

double g_K1 (const double *y)
{
    const double V_old_ = y[0];

    return (((1.2*exp((((-V_old_)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((V_old_+9.0e+01)/6.0e+01)))));
}

double i_k (const double *y)
{
    const double V_old_ = y[0];

    return (g_K1(y) + g_K2(y))*(V_old_ + 100.000);
}

double i_leak (const double *y)
{
    const double V_old_ = y[0];

    return g_L*(V_old_ - E_L);
}

double dvdt (double *y_star, double stim_current)
{
    double calc_I_stim = stim_current;

    return ( (-(i_na(y_star) + i_k(y_star) + i_leak(y_star)) + calc_I_stim ) ) / Cm;
}
// =================================================================================================

double beta_m (const double *y_old)
{
    const double V_old_ = y_old[0];

    return (((1.2e-01*(V_old_+8.0e+00))/(exp(((V_old_+8.0e+00)/5.0e+00))-1.0e+00)));
}

double alpha_m (const double *y_old)
{
    const double V_old_ = y_old[0];

    return (((1.0e-01*((-V_old_)-4.8e+01))/(exp((((-V_old_)-4.8e+01)/1.5e+01))-1.0e+00)));
}

double dmdt (double *y_old)
{
    const double m_old_ = y_old[1];

    return (alpha_m(y_old) *(1.00000 - m_old_) -  (beta_m(y_old)*m_old_) );
}
// =================================================================================================

double beta_h (const double *y_old)
{
    const double V_old_ = y_old[0];

    return ((1.0/(1.0e+00+exp((((-V_old_)-4.2e+01)/1.0e+01)))));
}

double alpha_h (const double *y_old)
{
    const double V_old_ = y_old[0];

    return ((1.7e-01*exp((((-V_old_)-9.0e+01)/2.0e+01))));
}

double dhdt (double *y_old)
{
    const double h_old_ = y_old[2];

    return (alpha_h(y_old) * (1.00000 - h_old_) -  (beta_h(y_old) * h_old_) );
}
// =================================================================================================

double beta_n (const double *y_old)
{
    const double V_old_ = y_old[0];

    return ((2.0e-03*exp((((-V_old_)-9.0e+01)/8.0e+01))));
}

double alpha_n (const double *y_old)
{
    const double V_old_ = y_old[0];

    return (((1.0e-04*((-V_old_)-5.0e+01))/(exp((((-V_old_)-5.0e+01)/1.0e+01))-1.0e+00)));
}

double dndt (double *y_old)
{
    const double n_old_ = y_old[3];

    return (alpha_n(y_old) * (1.00000 - n_old_) -  (beta_n(y_old) * n_old_) );
}
// =================================================================================================