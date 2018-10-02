#ifndef _CONFIG_PARSER_H_
#define _CONFIG_PARSER_H_

#include <cstdlib>
#include <cstdio>
#include <cstdbool>
#include <cctype>
#include <cmath>
#include <getopt.h>
#include "../../grid/grid.h"

#define SIGMA 1400              // Dummy values
#define STIM_OPT 1500           // Dummy values
#define BETA 4000               // Dummy values
#define CM 5000                 // Dummy values

struct user_options 
{
    int num_threads;
    bool num_threads_was_set;
    double final_time;
    bool final_time_was_set;
    double dt;
    bool dt_was_set;
    int print_rate;
    bool print_rate_was_set;
    char *out_dir_name;
    bool out_dir_name_was_set;
    int steady_state_print_rate;
    bool steady_state_print_rate_was_set;
    char *out_steady_state_dir;
    bool out_steady_state_dir_was_set;
    char *out_steady_state_filename;
    bool out_steady_state_filename_was_set;
    char *in_steady_state_filename;
    bool in_steady_state_filename_was_set;
    char *model_file_path;      
    bool model_file_path_was_set;
    bool gpu;                      
    bool gpu_was_set;
    bool use_steady_state;
    bool use_steady_state_was_set;
    int gpu_id;                    
    bool gpu_id_was_set;
    double sigma;
    bool sigma_was_set;
    double beta;
    bool beta_was_set;
    double cm;
    bool cm_was_set;

    char *config_file;

    struct stim_config_hash *stim_configs;
    struct purkinje_config *purkinje_config;
    
};

struct user_options * new_user_options();
void get_config_file (int argc, char **argv, struct user_options *user_args);
int parse_config_file (void *user, const char *section, const char *name, const char *value);
void parse_options (int argc, char **argv, struct user_options *user_args);
void display_usage (char **argv);
void issue_overwrite_warning (const char *var, const char *old_value, const char *new_value, const char *config_file) ;

void configure_grid_from_options(struct grid* grid, struct user_options *options);

#endif