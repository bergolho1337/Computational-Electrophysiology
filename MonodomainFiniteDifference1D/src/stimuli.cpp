#include "../include/stimuli.h"

struct stim_config* new_stim_config ()
{
	struct stim_config *stim = (struct stim_config*)malloc(sizeof(struct stim_config));

	return stim;
}

void free_stim_config (struct stim_config *stim)
{
	free(stim);
}

void configure_stimulus_from_options (struct stim_config *stim, struct user_options *options)
{
	stim->stim_start = options->stim_start;
	stim->stim_duration = options->stim_duration;
	stim->start_period = options->start_period;
	stim->end_period = options->end_period;
	stim->period_step = options->period_step;
	stim->n_cycles = options->n_cycles;

	print_stim_config(stim);
}

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

void print_stim_config (struct stim_config *stim)
{
	printf("%s\n",PRINT_LINE);
	printf("[Stimulus] stim_start = %.10lf\n",stim->stim_start);
	printf("[Stimulus] stim_duration = %.10lf\n",stim->stim_duration);
	printf("[Stimulus] start_period = %.10lf\n",stim->start_period);
	printf("[Stimulus] end_period = %.10lf\n",stim->end_period);
	printf("[Stimulus] period_step = %.10lf\n",stim->period_step);
	printf("[Stimulus] n_cycles = %d\n",stim->n_cycles);
	printf("%s\n",PRINT_LINE);
}
