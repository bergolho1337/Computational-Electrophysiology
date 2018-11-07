#include "fitzhugh_1961.h"

extern "C" GET_CELL_MODEL_DATA (get_cell_model_data)
{
    cell_model->initial_v = INITIAL_V;
    cell_model->number_of_ode_equations = NEQ;
}

extern "C" SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) 
{
    for (int i = 0; i < num_volumes; i++)
    {
        volumes[i].y_old = (double*)calloc(NEQ,sizeof(double));
        volumes[i].y_star = (double*)calloc(NEQ,sizeof(double));
        volumes[i].y_new = (double*)calloc(NEQ,sizeof(double));
        
        volumes[i].y_old[0] = 0.0;        // V
        volumes[i].y_old[1] = 0.0;          // h

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
    const double h_old_ = y_old[1];
    const double V_star_ = y_star[0];
    const double h_star_ = y_star[1];

    //Parameters
    const double alpha = -0.100000000000000e+00f;
    const double gamma = 3.000000000000000e+00f;
    const double epsilon = 5.000000000000000e-03f;

    double calc_I_stim = stim_current;

    rDY_[0] = (( V_star_*(V_star_ - alpha)*(1.00000 - V_star_) - h_star_) + calc_I_stim);
    rDY_[1] = epsilon*(V_old_ -  gamma*h_old_);

}
