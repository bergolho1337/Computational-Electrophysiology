#include "../include/stimulus.h"

Stimulus::Stimulus (Options *user_options)
{
    stim_current = user_options->stim_current;
    stim_start = user_options->stim_start;
    stim_duration = user_options->stim_duration;
    start_period = user_options->start_period;
    end_period = user_options->end_period;
    period_step = user_options->period_step;
    n_cycles = user_options->n_cycles;
    id_limit = user_options->id_limit;
}

void Stimulus::print ()
{
    cout << "/////////////////////////////////////////" << endl;
    cout << "Stim current = " << stim_current << endl;
    cout << "Stim start = " << stim_start << endl;
    cout << "Stim duration = " << stim_duration << endl;
    cout << "Start period = " << start_period << endl;
    cout << "End period = " << end_period << endl;
    cout << "Period step = " << period_step << endl;
    cout << "Number of cycles = " << n_cycles << endl;
    cout << "Id limit = " << id_limit << endl;
    cout << "/////////////////////////////////////////" << endl;
}