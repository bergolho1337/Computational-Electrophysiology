#include "stimuli.h"

struct stim_config* new_stim_config ()
{
    struct stim_config *config = (struct stim_config*)malloc(sizeof(struct stim_config));

    config->spatial_stim_currents = NULL;

    config->stim_current_was_set = false;
    config->stim_duration_was_set = false;
    config->stim_start_was_set = false;
    config->n_cycles_was_set = false;
    config->period_step_was_set = false;
    config->start_period_was_set = false;
    config->end_period_was_set = false;
    config->id_limit_was_set = false;

    return config;
}