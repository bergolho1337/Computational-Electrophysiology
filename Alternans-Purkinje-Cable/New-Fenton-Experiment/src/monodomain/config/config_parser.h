#ifndef FENTON_CONFIG_PARSER_H
#define FENTON_CONFIG_PARSER_H

#include <getopt.h>

#include <iostream>
#include <string>

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cstdbool>
#include <cctype>
#include <cmath>
//#include "../../alg/grid/grid.h"

#define SIGMA_X 1400
#define STIM_OPT 2000
#define BETA 4000
#define CM 5000

struct user_options 
{
    // Main section
    uint32_t num_threads;                
    double dt;					    
    double final_time;				
    bool use_steady_state;          
    uint32_t print_rate;	            	
    uint32_t sst_rate;	            	
    char* network_filename;
    char* sst_filename;
    char* plot_filename;

    // Cell section
    double start_h;
    uint32_t num_div_cell;
    double start_diameter;
    double sigma_c;
    double G_gap;
    char* model_file_path;

    // Stimulus
    struct stim_config_hash *stim_configs;

    // Configuration file
    std::string config_file;

};


/* Display program usage, and exit.
 */
void display_usage ( char *argv[] );
struct user_options* new_user_options();
void get_config_file (int argc, char **argv, struct user_options *user_args);
int parse_config_file(void* user, const char* section, const char* name, const char* value);
void print_user_options (struct user_options *user_args);
/*

void configure_grid_from_options(struct grid* grid, struct user_options *options);
void free_user_options(struct user_options *s);

void issue_overwrite_warning (const char *var, const char *old_value, const char *new_value, const char *config_file);
*/
#endif /* FENTON_CONFIG_PARSER_H */
