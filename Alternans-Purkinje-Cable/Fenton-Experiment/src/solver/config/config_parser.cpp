#include "config_parser.h"

void print_user_options (struct user_options *options)
{
    cout << PRINT_LINE << endl;
    cout << "num_threads = " << options->num_threads << endl;
    cout << "final_time = " << options->final_time << endl;
    cout << "dt = " << options->dt << endl;
    cout << "output_directory = " << options->out_dir_name << endl;
    cout << "print_rate = " << options->print_rate << endl;
    cout << "output_steady_state_directory = " << options->out_steady_state_dir << endl;
    cout << "steady_state_filename = " << options->steady_state_filename << endl;
    cout << "sigma = " << options->sigma << endl;
    cout << "cm = " << options->cm << endl;
    cout << "beta = " << options->beta << endl;
    cout << "config_file = " << options->config_file << endl;
    cout << PRINT_LINE << endl;
    cout << "name = " << options->purkinje_config->name << endl;
    cout << "network_file = " << options->purkinje_config->network_filename << endl;
    cout << "start_discretization = " << options->purkinje_config->start_h << endl;
    cout << PRINT_LINE << endl;
    cout << "stim_start = " << options->stim_configs->stim_start << endl;
    cout << "stim_current = " << options->stim_configs->stim_current << endl;
    cout << "stim_duration = " << options->stim_configs->stim_duration << endl;
    cout << "period_step = " << options->stim_configs->period_step << endl;
    cout << "start_period = " << options->stim_configs->start_period << endl;
    cout << "end_period = " << options->stim_configs->end_period << endl;
    cout << "n_cycles = " << options->stim_configs->n_cycles << endl;
    cout << "id_limit = " << options->stim_configs->id_limit << endl;
    cout << PRINT_LINE << endl;
}

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *options = (struct user_options*)malloc(sizeof(struct user_options));

    options->num_threads = 1;
    options->num_threads_was_set = false;

    options->out_dir_name = NULL;
    options->out_dir_name_was_set = false;

    options->final_time = 10.0;
    options->final_time_was_set = false;

    options->print_rate = 1;
    options->print_rate_was_set = false;

    options->out_steady_state_dir = NULL;
    options->out_steady_state_dir_was_set = false;

    options->dt = 0.01;
    options->dt_was_set = false;

    options->steady_state_filename = NULL;
    options->steady_state_filename_was_set = false;

    options->sigma = 0.0000176;
    options->sigma_was_set = false;

    options->cm = 1.2;
    options->cm_was_set = false;

    options->beta = 0.14;
    options->beta_was_set = false;

    options->config_file = NULL;

    options->purkinje_config = NULL;
    options->stim_configs = NULL;

    return options;
}