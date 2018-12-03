#include "../include/stimuli.h"

double get_spatial_stim_currents (const double x)
{
	static const double xmin = 0.0;
	static const double xmax = 0.2;

	static const double stim_current = 200.0;

	if (x >= xmin && x <= xmax)
		return stim_current;
	else
		return 0.0;
		
}

void compute_stimulus (double *stims, const double cur_time, const int np, const double dx)
{
	// Stimulus for the Noble 1962 celular model
	static const double stim_start = 0.0;
	static const double stim_duration = 5.0;
	static const double start_period = 300.0;
	static const double end_period = 300.0;
	static const double period_step = 100.0;
	static const int n_cycles = 10;

	double time = cur_time;
	double new_time, stim_period;
	
	new_time = 0.0;
        // New Jhonny stimulus protocol for alternans simulations ...
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
            //#pragma omp parallel for
            for (int i = 0; i < np; i++) 
            {
		double x = i*dx;
                stims[i] = get_spatial_stim_currents(x);
            }
        }
	else
	{
	    //#pragma omp parallel for
	    for (int i = 0; i < np; i++)
		stims[i] = 0.0;
	}
        time = cur_time;
}


