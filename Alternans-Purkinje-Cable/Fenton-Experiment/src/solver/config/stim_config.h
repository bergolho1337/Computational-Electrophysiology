#include <cstdio>
#include <cstdlib>
#include <cstring>

struct stim_config
{
    double stim_start;
    bool stim_start_was_set;
    double stim_duration;
    bool stim_duration_was_set;
    double stim_current;
    bool stim_current_was_set;
    // Variables related to Jhonny's stimulus protocol ...
    int n_cycles;
    bool n_cycles_was_set;
    double start_period;
    bool start_period_was_set;
    double end_period;
    bool end_period_was_set;
    double period_step;
    bool period_step_was_set;
    // Berg's variables
    int id_limit;
    bool id_limit_was_set;

    double *spatial_stim_currents;

};

struct stim_config* new_stim_config ();