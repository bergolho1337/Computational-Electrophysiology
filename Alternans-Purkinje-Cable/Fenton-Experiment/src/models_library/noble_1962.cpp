#include "noble_1962.h"

extern "C" GET_CELL_MODEL_DATA (get_cell_model_data)
{
    cell_model->initial_v = INITIAL_V;
    cell_model->number_of_ode_equations = NEQ;

    for (int i = 0; i < num_volumes; i++)
    {
        volumes[i].y_old = (double*)calloc(NEQ,sizeof(double));
        volumes[i].y_star = (double*)calloc(NEQ,sizeof(double));
        volumes[i].y_new = (double*)calloc(NEQ,sizeof(double));
    }
}

extern "C" SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) 
{
    for (int i = 0; i < num_volumes; i++)
    {
        
        // Old initial conditions
        //volumes[i].y_old[0] = -75.5344986658;   // V
        //volumes[i].y_old[1] = 0.0605467272;     // m
        //volumes[i].y_old[2] = 0.7259001355;     // h       
        //volumes[i].y_old[3] = 0.4709239708;     // n

        // New initial conditions
        volumes[i].y_old[0] = -87.0;   // V
        volumes[i].y_old[1] = 0.01;     // m
        volumes[i].y_old[2] = 0.8;     // h       
        volumes[i].y_old[3] = 0.01;     // n

    }
}

extern "C" SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) 
{
    int i;

    //#pragma omp parallel for
    for (i = 0; i < num_volumes; i++)
    {
        solve_model_ode_cpu(dt,volumes[i],stim_currents[i]);
    }
    
}

void solve_model_ode_cpu(double dt, struct control_volume &volume,\
                         double stim_current)  
{

    double *y_old = volume.y_old;
    double *y_star = volume.y_star;
    double *y_new = volume.y_new;

    double rDY[NEQ];
    RHS_cpu(rDY,y_old,y_star,stim_current);

    for (int i = 0; i < NEQ; i++)
        y_new[i] = dt*rDY[i] + y_star[i];

}

void RHS_cpu(double *rDY_, const double *y_old, const double *y_star, double stim_current) 
{

    //State variables
    const double V_old_ = y_old[0];
    const double m_old_ = y_old[1];
    const double h_old_ = y_old[2];
    const double n_old_ = y_old[3];

    const double V_star_ = y_star[0];
    const double m_star_ = y_star[1];
    const double h_star_ = y_star[2];
    const double n_star_ = y_star[3];

    //Parameters
    const double g_Na_max = 4.0e+02;
    const double E_Na = 4.0e+01;
    const double g_L = 7.5e-02;
    const double E_L = -6.0e+01;
    const double CM = 1.2e+01;

    // Algebraics
    double gna = ((pow(m_star_,3.0e+00)*h_star_*g_Na_max));
    double ina = (gna+1.4e-01)*(V_star_-E_Na);
    double gk1 = (((1.2*exp((((-V_star_)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((V_star_+9.0e+01)/6.0e+01)))));
    double gk2 = ((1.2*pow(n_star_,4.0e+00)));
    double ik = (((gk1+gk2)*(V_star_+1.0e+02)));
    double ileak = ((g_L*(V_star_-E_L)));
    double istim = stim_current;
    
    double beta_m = (((1.2e-01*(V_old_+8.0e+00))/(exp(((V_old_+8.0e+00)/5.0e+00))-1.0e+00)));
    double alpha_m = (((1.0e-01*((-V_old_)-4.8e+01))/(exp((((-V_old_)-4.8e+01)/1.5e+01))-1.0e+00)));
    double beta_h = ((1.0/(1.0e+00+exp((((-V_old_)-4.2e+01)/1.0e+01)))));
    double alpha_h = ((1.7e-01*exp((((-V_old_)-9.0e+01)/2.0e+01))));
    double beta_n = ((2.0e-03*exp((((-V_old_)-9.0e+01)/8.0e+01))));
    double alpha_n = (((1.0e-04*((-V_old_)-5.0e+01))/(exp((((-V_old_)-5.0e+01)/1.0e+01))-1.0e+00)));


    rDY_[0] = ((-(ina+ik+ileak)+istim)/CM);
    rDY_[1] = ((alpha_m*(1.0e+00-m_old_))-(beta_m*m_old_));
    rDY_[2] = ((alpha_h*(1.0e+00-h_old_))-(beta_h*h_old_));
    rDY_[3] = (alpha_n*(1.0e+00-n_old_))-(beta_n*n_old_);

}
