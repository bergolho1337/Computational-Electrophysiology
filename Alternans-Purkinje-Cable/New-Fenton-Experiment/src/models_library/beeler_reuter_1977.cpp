#include "beeler_reuter_1977.h"

extern "C" GET_CELL_MODEL_DATA(get_cell_model_data) 
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
        
        volumes[i].y_old[0] = -84.624;        // V
        volumes[i].y_old[1] = 0.011;          // m
        volumes[i].y_old[2] = 0.988;          // h
        volumes[i].y_old[3] = 0.975;          // j
        volumes[i].y_old[4] = 1e-4;           // Cai
        volumes[i].y_old[5] = 0.003;          // d
        volumes[i].y_old[6] = 0.994;          // f
        volumes[i].y_old[7] = 0.0001;         // x1

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

    double rDY[NEQ];
    RHS_cpu(dt,rDY,volume.y_old,volume.y_star,stim_current);

    // Old Euler code
    for (int i = 0; i < NEQ; i++)
    {
        volume.y_new[i] = dt*rDY[i] + volume.y_star[i];
    }
        

    /*
    // Forward Euler
    y_new[0] = dt*rDY[0] + y_old[0];
    y_new[4] = dt*rDY[4] + y_old[4];

    // Rush-Larsen
    y_new[1] = rDY[1];
    y_new[2] = rDY[2];
    y_new[3] = rDY[3];
    y_new[5] = rDY[5];
    y_new[6] = rDY[6];
    y_new[7] = rDY[7];
    */

}

void RHS_cpu(const double dt,double *rDY_, const double *y_old, const double *y_star, double stim_current) 
{

    // State variables
    const double V_old_ = y_old[0];
    const double m_old_ = y_old[1];
    const double h_old_ = y_old[2];
    const double j_old_ = y_old[3];
    const double Cai_old_ = y_old[4];
    const double d_old_ = y_old[5];
    const double f_old_ = y_old[6];
    const double x1_old_ = y_old[7];

    const double V_star_ = y_star[0];
    const double m_star_ = y_star[1];
    const double h_star_ = y_star[2];
    const double j_star_ = y_star[3];
    const double Cai_star_ = y_star[4];
    const double d_star_ = y_star[5];
    const double f_star_ = y_star[6];
    const double x1_star_ = y_star[7];

    // Constants
    const double C = 0.01;
    const double g_na = 4e-2;
    const double E_na = 50;
    const double g_nac = 3e-5;
    const double g_s = 9e-4;

    // Algebraics
    double alpha_m = ( - 1.00000*(V_old_+47.0000))/(exp( - 0.100000*(V_old_+47.0000)) - 1.00000);
    double beta_m =  40.0000*exp( - 0.0560000*(V_old_+72.0000));
    double alpha_h =  0.126000*exp( - 0.250000*(V_old_+77.0000));
    double beta_h = 1.70000/(exp( - 0.0820000*(V_old_+22.5000))+1.00000);
    double alpha_j = ( 0.0550000*exp( - 0.250000*(V_old_+78.0000)))/(exp( - 0.200000*(V_old_+78.0000))+1.00000);
    double beta_j = 0.300000/(exp( - 0.100000*(V_old_+32.0000))+1.00000);
    double alpha_d = ( 0.0950000*exp(- (V_old_ - 5.00000)/100.000))/(1.00000+exp(- (V_old_ - 5.00000)/13.8900));
    double beta_d = ( 0.0700000*exp(- (V_old_+44.0000)/59.0000))/(1.00000+exp((V_old_+44.0000)/20.0000));
    double alpha_f = ( 0.0120000*exp(- (V_old_+28.0000)/125.000))/(1.00000+exp((V_old_+28.0000)/6.67000));
    double beta_f = ( 0.00650000*exp(- (V_old_+30.0000)/50.0000))/(1.00000+exp(- (V_old_+30.0000)/5.00000));
    double alpha_x1 = ( 0.000500000*exp((V_old_+50.0000)/12.1000))/(1.00000+exp((V_old_+50.0000)/17.5000));
    double beta_x1 = ( 0.00130000*exp(- (V_old_+20.0000)/16.6700))/(1.00000+exp(- (V_old_+20.0000)/25.0000));
    double E_s = - 82.3000 -  13.0287*log( Cai_old_*0.00100000);
    double i_s =  g_s*d_old_*f_old_*(V_old_ - E_s);

    double i_na =  ( g_na*pow(m_star_, 3.00000)*h_star_*j_star_+g_nac)*(V_star_ - E_na);
    double i_x1 = ( x1_star_*0.00800000*(exp( 0.0400000*(V_star_+77.0000)) - 1.00000))/exp( 0.0400000*(V_star_+35.0000));
    double i_k1 =  0.00350000*(( 4.00000*(exp( 0.0400000*(V_star_+85.0000)) - 1.00000))/(exp( 0.0800000*(V_star_+53.0000))+exp( 0.0400000*(V_star_+53.0000)))+( 0.200000*(V_star_+23.0000))/(1.00000 - exp( - 0.0400000*(V_star_+23.0000))));
    double i_stim = stim_current;

    // Rates
    rDY_[0] = (i_stim - (i_na+i_s+i_x1+i_k1))/C;
    rDY_[1] = alpha_m*(1.00000 - m_old_) -  beta_m*m_old_;
    rDY_[2] = alpha_h*(1.00000 - h_old_) -  beta_h*h_old_;
    rDY_[3] = alpha_j*(1.00000 - j_old_) -  beta_j*j_old_;
    rDY_[4] = ( - 0.0100000*i_s)/1.00000+ 0.0700000*(0.000100000 - Cai_old_);
    rDY_[5] = alpha_d*(1.00000 - d_old_) -  beta_d*d_old_;
    rDY_[6] = alpha_f*(1.00000 - f_old_) -  beta_f*f_old_;
    rDY_[7] = alpha_x1*(1.00000 - x1_old_) -  beta_x1*x1_old_;

}