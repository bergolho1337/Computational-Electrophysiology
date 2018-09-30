#include "config_parser.h"

struct user_options * new_user_options()
{
    struct user_options *user_args = (struct user_options*)malloc(sizeof(struct user_options));

    user_args->out_dir_name = NULL;
    user_args->out_dir_name_was_set = false;

    user_args->out_steady_state_dir = NULL;
    user_args->out_steady_state_dir_was_set = false;

    user_args->steady_state_filename = NULL;
    user_args->steady_state_filename_was_set = false;

    user_args->num_threads = 1;
    user_args->num_threads_was_set = false;

    user_args->final_time = 10.0;
    user_args->final_time_was_set = false;

    user_args->print_rate = 1;
    user_args->print_rate_was_set = false;

    user_args->model_file_path = NULL;
    user_args->model_file_path_was_set = false;

    user_args->gpu = false;
    user_args->gpu_was_set = false;

    user_args->gpu_id = 0;
    user_args->gpu_id_was_set = false;

    user_args->dt = 0.01;
    user_args->dt_was_set = false;

    user_args->config_file = NULL;

    user_args->sigma = 0.0000176;
    user_args->sigma_was_set = false;

    user_args->beta = 0.14;
    user_args->beta_was_set = false;

    user_args->cm = 1.2;
    user_args->cm_was_set = false;

    return user_args;
}