#include "fitzhugh_1961.h"

GET_CELL_MODEL_DATA(init_cell_model_data) 
{

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) 
{

    sv[0] = 0.000000f; //V millivolt 
    sv[1] = 0.000000f; //h dimensionless 
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) 
{

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

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  
{

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current);

    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current) 
{

    //State variables
    const real V_old_ = sv[0];
    const real h_old_ = sv[1];

    //Parameters
    const real alpha = -0.100000000000000e+00f;
    const real gamma = 3.000000000000000e+00f;
    const real epsilon = 5.000000000000000e-03f;

    real calc_I_stim = stim_current;

    rDY_[0] = (( V_old_*(V_old_ - alpha)*(1.00000 - V_old_) - h_old_) + calc_I_stim);
    rDY_[1] = epsilon*(V_old_ -  gamma*h_old_);

}
