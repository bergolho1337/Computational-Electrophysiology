#include <cstring>
#include <cassert>

#include "config_parser.h"
#include "../../ini_parser/ini_file_sections.h"
#include "../../utils/logfile_utils.h"

#include "stim_config_hash.h"
#include "purkinje_config.h"

static const char *opt_string =   "c:o:n:p:t:e:f:k:g:G:h";

static const struct option long_options[] = {
        { "config_file", required_argument, NULL, 'c' },
        { "output_dir", required_argument , NULL, 'o' },
        { "num_threads", required_argument, NULL, 'n' },
        { "print_rate",required_argument , NULL, 'p' },
        { "dt", required_argument, NULL, 'e' },
        { "simulation_time", required_argument, NULL, 'f' },
        { "model_file_path", required_argument, NULL, 'k'},
        { "gpu_id", required_argument, NULL, 'G'},
        { "use_gpu", required_argument, NULL, 'g' },
        { "sigma", required_argument, NULL, SIGMA},
        { "beta", required_argument, NULL, BETA},
        { "cm", required_argument, NULL, CM},
        { "stimulus", required_argument, NULL, STIM_OPT},
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
};

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

void get_config_file (int argc, char **argv, struct user_options *user_args) 
{
    int opt = 0;

    int option_index;

    opt = getopt_long (argc, argv, opt_string, long_options, &option_index);

    while (opt != -1) 
    {
        switch (opt) 
        {
            case 'c':
                user_args->config_file = optarg;
                return;
            default:
                break;
        }
        opt = getopt_long (argc, argv, opt_string, long_options, &option_index);
    }

    // We reset the index after parsing the config_file
    optind = 1;
}

int parse_config_file (void *user, const char *section, const char *name, const char *value) 
{
    struct user_options *pconfig = (struct user_options *) user;

    if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "num_threads")) 
    {
        pconfig->num_threads = (int)strtol (value, NULL, 10);
        pconfig->num_threads_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "dt")) 
    {
        pconfig->dt = strtod(value, NULL);
        pconfig->dt_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "simulation_time")) 
    {
        pconfig->final_time = strtod(value, NULL);
        pconfig->final_time_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "sigma")) 
    {
        pconfig->sigma = strtod(value, NULL);
        pconfig->sigma_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "beta")) 
    {
        pconfig->beta = strtod(value, NULL);
        pconfig->beta_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "cm")) 
    {
        pconfig->cm = strtod(value, NULL);
        pconfig->cm_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "print_rate")) 
    {
        pconfig->print_rate = (int)strtol (value, NULL, 10);
        pconfig->print_rate_was_set = true;
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "output_dir")) 
    {
        pconfig->out_dir_name = strdup(value);
        pconfig->out_dir_name_was_set = true;
    }
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "output_steady_state_dir")) 
    {
        pconfig->out_steady_state_dir = strdup(value);
        pconfig->out_steady_state_dir_was_set = true;
    }
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "steady_state_filename")) 
    {
        pconfig->steady_state_filename = strdup(value);
        pconfig->steady_state_filename_was_set = true;
    } 
    else if (MATCH_SECTION(ODE_SECTION)) 
    {
        if(MATCH_NAME("use_gpu")) 
        {
            if (strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) 
            {
                pconfig->gpu = true;
            } 
            else 
            {
                pconfig->gpu = false;
            }
            pconfig->gpu_was_set = true;
        } 
        else if (MATCH_NAME ("gpu_id")) 
        {
            pconfig->gpu_id = (int)strtol (value, NULL, 10);
            pconfig->gpu_id_was_set = true;
        } 
        else if (MATCH_NAME ("library_file")) {
            pconfig->model_file_path = strdup(value);
            pconfig->model_file_path_was_set = true;
        }
    }
    else if(SECTION_STARTS_WITH(STIM_SECTION)) 
    {

        if(pconfig->stim_configs == NULL) 
        {
            pconfig->stim_configs = stim_config_hash_create();
        }

        struct stim_config *tmp = stim_config_hash_search(pconfig->stim_configs, section);

        if( tmp == NULL) 
        {
            tmp = new_stim_config();
            stim_config_hash_insert(pconfig->stim_configs, section, tmp);
        }

        if(MATCH_NAME("stim_start")) 
        {
            tmp->stim_start = (float)strtod(value, NULL);
            tmp->stim_start_was_set = true;
        } 
        else if(MATCH_NAME("stim_duration")) 
        {
            tmp->stim_duration = (float)strtod(value, NULL);
            tmp->stim_duration_was_set = true;
        }
        else if(MATCH_NAME("stim_current")) 
        {
            tmp->stim_current = (float)strtod(value, NULL);
            tmp->stim_current_was_set = true;
        }
        else if(MATCH_NAME("n_cycles")) 
        {
            tmp->n_cycles = (int)strtol(value, NULL, 10);
            tmp->n_cycles_was_set = true;
        }
        else if(MATCH_NAME("start_period")) 
        {
            tmp->start_period = (float)strtod(value, NULL);
            tmp->start_period_was_set = true;
        }
        else if(MATCH_NAME("end_period")) 
        {
            tmp->end_period = (float)strtod(value, NULL);
            tmp->end_period_was_set = true;
        }   
        else if(MATCH_NAME("period_step")) 
        {
            tmp->period_step = (float)strtod(value, NULL);
            tmp->period_step_was_set = true;
        }
        else if(MATCH_NAME("function")) 
        {
            tmp->config_data.function_name = strdup(value);
            tmp->config_data.function_name_was_set = true;
        } 
        else if(MATCH_NAME("library_file")) 
        {
            tmp->config_data.library_file_path = strdup(value);
            tmp->config_data.library_file_path_was_set = true;
        }
        else if(MATCH_NAME("id_limit")) 
        {
            tmp->id_limit = (int)strtol(value, NULL, 10);
            tmp->config_data.library_file_path_was_set = true;
        }
    }
    else if(MATCH_SECTION(PURKINJE_SECTION)) 
    {

        if(pconfig->purkinje_config == NULL) 
        {
            pconfig->purkinje_config = new_purkinje_config();
        }

        if (MATCH_NAME ("start_discretization")) 
        {
            pconfig->purkinje_config->start_h = strtod (value, NULL);
            pconfig->purkinje_config->start_h_was_set = true;
        }
        else if (MATCH_NAME ("diameter")) 
        {
            pconfig->purkinje_config->diameter = strtod (value, NULL);
            pconfig->purkinje_config->diameter_was_set = true;
        }
        else if(MATCH_NAME("name")) 
        {
            pconfig->purkinje_config->network_name = strdup(value);
            pconfig->purkinje_config->network_name_was_set = true;
        }
        else if(MATCH_NAME("network_file")) 
        {
            pconfig->purkinje_config->network_filename = strdup(value);
            pconfig->purkinje_config->network_filename_was_set = true;
        }
        else if(MATCH_NAME("function")) 
        {
            pconfig->purkinje_config->config_data.function_name = strdup(value);
            pconfig->purkinje_config->config_data.function_name_was_set = true;
        }
        else if(MATCH_NAME("library_file")) 
        {
            pconfig->purkinje_config->config_data.library_file_path = strdup(value);
            pconfig->purkinje_config->config_data.library_file_path_was_set = true;
        }
    }
    else 
    {
        fprintf(stderr, "Invalid name %s in section %s on the config file!\n", name, section);
        return 0;
    }
    return 1;
}

void parse_options (int argc, char **argv, struct user_options *user_args) 
{

    int opt = 0;
    int option_index;

    opt =  getopt_long_only (argc, argv, opt_string, long_options, &option_index);
    char old_value[32];

    while (opt != -1) 
    {
        switch (opt) 
        {
    
            case 'p':
                if (user_args->print_rate_was_set) 
                {
                    sprintf (old_value, "%d", user_args->print_rate);
                    issue_overwrite_warning ("print_rate", old_value, optarg, user_args->config_file);
                }
                user_args->print_rate = (int)strtol (optarg, NULL, 10);
                break;
            case 'o':
                if (user_args->out_dir_name_was_set) 
                {
                    if (user_args->out_dir_name) 
                    {
                        issue_overwrite_warning ("output_dir", user_args->out_dir_name, optarg, user_args->config_file);
                    } 
                    else 
                    {
                        issue_overwrite_warning ("output_dir", "No Save", optarg, user_args->config_file);
                    }
                }
                free(user_args->out_dir_name);
                user_args->out_dir_name = strdup(optarg);

                break;
            case 'k':
                if (user_args->model_file_path_was_set) 
                {
                    if (user_args->model_file_path) 
                    {
                        issue_overwrite_warning ("model_file_path", user_args->model_file_path, optarg,
                                                 user_args->config_file);
                    } 
                    else 
                    {
                        issue_overwrite_warning ("model_file_path", "No Save", optarg, user_args->config_file);
                    }
                }
                free(user_args->model_file_path);
                user_args->model_file_path = strdup(optarg);

                break;
            case 'f':
                if (user_args->final_time_was_set) 
                {
                    sprintf (old_value, "%lf", user_args->final_time);
                    issue_overwrite_warning ("simulation_time", old_value, optarg, user_args->config_file);
                }
                user_args->final_time = strtod (optarg, NULL);

                break;
            case 'n':
                if (((int)strtol (optarg, NULL, 10)) > 0) 
                {
                    if (user_args->num_threads_was_set) 
                    {
                        sprintf (old_value, "%d", user_args->num_threads);
                        issue_overwrite_warning ("nu"
                                                         "m_threads", old_value, optarg, user_args->config_file);
                    }
                    user_args->num_threads = (int)strtol (optarg, NULL, 10);
                }
                break;
            case 'z':
                if (user_args->dt_was_set) 
                {
                    sprintf (old_value, "%lf", user_args->dt);
                    issue_overwrite_warning ("dt", old_value, optarg, user_args->config_file);
                }
                user_args->dt = strtod (optarg, NULL);
                break;
            case SIGMA:
                if (user_args->sigma_was_set) 
                {
                    sprintf (old_value, "%lf", user_args->sigma);
                    issue_overwrite_warning ("sigma_x", old_value, optarg, user_args->config_file);
                }
                user_args->sigma = strtod (optarg, NULL);
                break;
            case BETA:
                if (user_args->beta_was_set) 
                {
                    sprintf (old_value, "%lf", user_args->beta);
                    issue_overwrite_warning ("beta", old_value, optarg, user_args->config_file);
                }
                user_args->beta = strtod (optarg, NULL);
                break;
            case CM:
                if (user_args->cm) 
                {
                    sprintf (old_value, "%lf", user_args->cm);
                    issue_overwrite_warning ("cm", old_value, optarg, user_args->config_file);
                }
                user_args->cm = strtod (optarg, NULL);
                break;
                /*
            case STIM_OPT:
                if(user_args->stim_configs == NULL) 
                {
                    print_to_stdout_and_file("Creating new stim config from command line!\n");
                    user_args->stim_configs = stim_config_hash_create();
                }
                set_stim_config(optarg, user_args->stim_configs, user_args->config_file );
                break;
                */
            case 'h': /* fall-through is intentional */
            case '?':
                display_usage(argv);
                break;
            default:
                /* You won't actually get here. */
                break;
        }

        opt = getopt_long (argc, argv, opt_string, long_options, &option_index);
    }
}

void issue_overwrite_warning (const char *var, const char *old_value, const char *new_value, const char *config_file) 
{
    fprintf (stderr,
             "WARNING: option %s was set in the file %s to %s and is being overwritten "
                     "by the command line flag to %s!\n",
             var, config_file, old_value, new_value);
}

void display_usage (char **argv) 
{

    //TODO: help for domain, extra_data and stimuls flags
    printf ("Usage: %s [options]\n\n", argv[0]);
    printf ("Options:\n");
    printf ("--config_file | -c [configuration_file_path]. Simulation configuration file. Default NULL.\n");
    printf ("--simulation_time | -f [simulation final time]. Simulation final time. Default 10.\n");
    printf ("--output_dir | -o [output-dir]. Simulation output directory. If not set, the simulation will not be saved "
                    "to disk. \n");
    printf ("--print_rate | -p [output-print-rate]. Output print rate (in number of iterations). Default: 1 \n");
    printf ("--dt | -z [dt]. Simulation time discretization (PDE && ODE). Default: 0.01 \n");
    printf ("--beta. Value of beta for simulation (Default: 0.14 \n");
    printf ("--cm. Value cm (Default: 1.2 \n");
    printf ("--num_threads | -n [num-threads]. Solve using OpenMP. Default: 1 \n");
    printf ("--model_file_path | -k [.so file path], Path of the .so representing the cell model and the ode solver. "
                    "Default: NULL \n");
    printf ("--help | -h. Shows this help and exit \n");
    exit (EXIT_FAILURE);
}

void configure_grid_from_options(struct grid* grid, struct user_options *options) 
{
    assert(grid);
    assert(options);
}