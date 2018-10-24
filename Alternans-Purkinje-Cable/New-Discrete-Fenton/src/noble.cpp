#include "../include/noble.h"

Stimulus *stim;

// New Jhonny stimulus protocol for alternans simulations ...
double I_Stim__Nob (int point, double t)
{
    double stim_current = stim->stim_current;
    double stim_start = stim->stim_start;
    double stim_duration = stim->stim_duration;
    double start_period = stim->start_period;
    double end_period = stim->end_period;
    double period_step = stim->period_step;
    int n_cycles = stim->n_cycles;
    int id_limit = stim->id_limit;

    double stim_period;
    double new_time = 0.0f;
    double time = t; 
    if (point <= id_limit)
    {
        
        for (double new_period = start_period; new_period >= end_period; new_period -= period_step)
        {
            if ( time >= new_time && (time < new_time + n_cycles*new_period || new_period == end_period) )
            {
                stim_period = new_period;
                time -= new_time;
                break;
            }
            new_time += n_cycles*new_period;

        }
        if( (time-floor(time/stim_period)*stim_period>=stim_start) && ( time - floor(time/stim_period)*stim_period <= stim_start + stim_duration ) )
        {
            return stim_current;
        }
        else
        {
            return 0.0f;
        }
    }
    else
    {
        return 0.0f;
    }
}

double I_Leak__Nob (double t, double vm, double m, double h, double n)
{
	return ((g_L*(vm-E_L)));
}

// Oscillations can be eliminated by adding 1.0e+02
double g_K2__Nob (double t, double vm, double m, double h, double n)
{
    return ((1.2*pow(n,4.0e+00)));
	//return ((1.2e+03*pow(y[3],4.0e+00))+1.0e-01);          // Without oscillations
}

// Change the first coefficient to 1.2e+03 to 1.3e+03 can also eliminate the oscillations
double g_K1__Nob (double t, double vm, double m, double h, double n)
{
	//return (((1.2*exp((((-vm)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((vm+9.0e+01)/6.0e+01)))));
    //return (((1.3*exp((((-vm)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((vm+9.0e+01)/6.0e+01)))));    // Original
    return (((1.2*exp((((-vm)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((vm+9.0e+01)/6.0e+01)))));
}

double I_K__Nob (double t, double vm, double m, double h, double n)
{
    return (((g_K1__Nob(t,vm,m,h,n)+g_K2__Nob(t,vm,m,h,n))*(vm+1.0e+02)));
}

double g_Na__Nob (double t, double vm, double m, double h, double n)
{
    return ((pow(m,3.0e+00)*h*g_Na_max));
}

// Change the coeffcient to 1.4e+02 to 1.225e+02 can also eliminate the oscillations
double I_Na__Nob (double t, double vm, double m, double h, double n)
{
    return ((g_Na__Nob(t,vm,m,h,n)+1.4e-01)*(vm-E_Na));
}

// Sodium current without oscillations
double I_Na_NoOscilation__Nob (double t, double vm, double m, double h, double n)
{
    return ((g_Na__Nob(t,vm,m,h,n)+1.220e-01)*(vm-E_Na));      // Without oscillations
}

double dvdt__Nob (int type, int point, double t, double y[])
{
    // Purkinje cell
    if (type == 0)
        return ((-(I_Na__Nob(t,y[0],y[1],y[2],y[3])+I_K__Nob(t,y[0],y[1],y[2],y[3])+I_Leak__Nob(t,y[0],y[1],y[2],y[3]))+I_Stim__Nob(point,t))/CM);
    // Myocardium cell
    else
        return ((-(I_Na_NoOscilation__Nob(t,y[0],y[1],y[2],y[3])+I_K__Nob(t,y[0],y[1],y[2],y[3])+I_Leak__Nob(t,y[0],y[1],y[2],y[3])))/CM);
    //return ((-(I_Na(t,y)+I_K(t,y)+I_Leak(t,y))/CM));
}

/* ====================================================================================================== */
// dm/dt

double beta_m__Nob (double t, double vm, double m, double h, double n)
{
    return (((1.2e-01*(vm+8.0e+00))/(exp(((vm+8.0e+00)/5.0e+00))-1.0e+00)));
}

double alpha_m__Nob (double t, double vm, double m, double h, double n)
{
	return (((1.0e-01*((-vm)-4.8e+01))/(exp((((-vm)-4.8e+01)/1.5e+01))-1.0e+00)));
}

double dmdt__Nob (int type, int point, double t, double y[])
{
	return ((alpha_m__Nob(t,y[0],y[1],y[2],y[3])*(1.0e+00-y[1]))-(beta_m__Nob(t,y[0],y[1],y[2],y[3])*y[1]));
}

/* ====================================================================================================== */
// dh/dt

double beta_h__Nob (double t, double vm, double m, double h, double n)
{
	return ((1.0/(1.0e+00+exp((((-vm)-4.2e+01)/1.0e+01)))));
}

double alpha_h__Nob (double t, double vm, double m, double h, double n)
{
	return ((1.7e-01*exp((((-vm)-9.0e+01)/2.0e+01))));
}

double dhdt__Nob (int type, int point, double t, double y[])
{
	return ((alpha_h__Nob(t,y[0],y[1],y[2],y[3])*(1.0e+00-y[2]))-(beta_h__Nob(t,y[0],y[1],y[2],y[3])*y[2]));
}

/* ====================================================================================================== */
// dn/dt

double beta_n__Nob (double t, double vm, double m, double h, double n)
{
	return ((2.0e-03*exp((((-vm)-9.0e+01)/8.0e+01))));
}

double alpha_n__Nob (double t, double vm, double m, double h, double n)
{
	return (((1.0e-04*((-vm)-5.0e+01))/(exp((((-vm)-5.0e+01)/1.0e+01))-1.0e+00)));
}

double dndt__Nob (int type, int point, double t, double y[])
{
	return ((alpha_n__Nob(t,y[0],y[1],y[2],y[3])*(1.0e+00-y[3]))-(beta_n__Nob(t,y[0],y[1],y[2],y[3])*y[3]));
}

/* ====================================================================================================== */

void set_celular_model (Stimulus *stim_config)
{
    stim = stim_config;
}