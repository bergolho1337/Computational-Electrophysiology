#include "config_parser.h"
#include "../../ini_parser/ini_file_sections.h"

#include "stim_config_hash.h"

static const struct option long_options[] = {
        { "config_file", required_argument, NULL, 'c' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
};

static const char *opt_string = "c:h";


struct user_options* new_user_options()
{
    struct user_options *user_args = (struct user_options *)malloc (sizeof (struct user_options));

    user_args->num_threads = 1;
    user_args->dt = 0.02;
    user_args->final_time = 10.0;
    user_args->use_steady_state = false;
    user_args->print_rate = 1;
    user_args->sst_rate = 1;

    user_args->stim_configs = NULL;

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
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "dt")) {
        pconfig->dt = strtod(value, NULL);
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "simulation_time")) 
    {
        pconfig->final_time = strtod(value, NULL);
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "use_steady_state")) 
    {
        if (strcmp(value, "true") == 0 || strcmp(value, "yes") == 0)
            pconfig->use_steady_state = true;
        else
            pconfig->use_steady_state = false;
    }
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "print_rate")) 
    {
        pconfig->print_rate = (int)strtol (value, NULL, 10);
    }
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "sst_rate")) 
    {
        pconfig->sst_rate = (int)strtol (value, NULL, 10);
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "network_filename")) 
    {
        pconfig->network_filename = strdup(value);
    } 
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "sst_filename")) 
    {
        pconfig->sst_filename = value;
    }
    else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "plot_filename")) 
    {
        pconfig->plot_filename = value;
    }

    else if (MATCH_SECTION_AND_NAME (CELL_SECTION, "start_h")) 
    {
        pconfig->start_h = strtod(value, NULL);
    }
    else if (MATCH_SECTION_AND_NAME (CELL_SECTION, "num_div_cell")) 
    {
        pconfig->num_div_cell = (int)strtol (value, NULL, 10);
    }
    else if (MATCH_SECTION_AND_NAME (CELL_SECTION, "start_diameter")) 
    {
        pconfig->start_diameter = strtod(value, NULL);
    } 
    else if (MATCH_SECTION_AND_NAME (CELL_SECTION, "sigma_c")) 
    {
        pconfig->sigma_c = strtod(value, NULL);
    }
    else if (MATCH_SECTION_AND_NAME (CELL_SECTION, "G_gap")) 
    {
        pconfig->G_gap = strtod(value, NULL);
    }
    else if (MATCH_SECTION_AND_NAME (CELL_SECTION, "library_file")) 
    {
        pconfig->model_file_path = value;
    }

    else if(SECTION_STARTS_WITH(STIM_SECTION)) 
    {

        if(pconfig->stim_configs == NULL) 
        {
            pconfig->stim_configs = stim_config_hash_create();
        }

        // Search if the stimulus protocol is already inserted in the hash table
        struct stim_config *tmp = stim_config_hash_search(pconfig->stim_configs, section);

        // Put into the new stimulus protocol configuration in the hash table
        if( tmp == NULL) 
        {
            tmp = new_stim_config();
            stim_config_hash_insert(pconfig->stim_configs, section, tmp);
        }

        if(MATCH_NAME("stim_start")) 
        {
            tmp->stim_start = (real)strtod(value, NULL);
            tmp->stim_start_was_set = true;
        } 
        else if(MATCH_NAME("stim_duration")) 
        {
            tmp->stim_duration = (real)strtod(value, NULL);
            tmp->stim_duration_was_set = true;
        }
        else if(MATCH_NAME("stim_current")) 
        {
            tmp->stim_current = (real)strtod(value, NULL);
            tmp->stim_current_was_set = true;
        }
        else if(MATCH_NAME("n_cycles")) 
        {
            tmp->n_cycles = (int)strtol(value, NULL, 10);
            tmp->n_cycles_was_set = true;
        }
        else if(MATCH_NAME("start_period")) 
        {
            tmp->start_period = (real)strtod(value, NULL);
            tmp->start_period_was_set = true;
        }
        else if(MATCH_NAME("end_period")) 
        {
            tmp->end_period = (real)strtod(value, NULL);
            tmp->end_period_was_set = true;
        }   
        else if(MATCH_NAME("period_step")) 
        {
            tmp->period_step = (real)strtod(value, NULL);
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
        else 
        {
            //name is a reserved word in stim config
            if(MATCH_NAME("name")) 
            {
                fprintf(stderr, "name is a reserved word and should not be used inside a stimulus config section. Found in %s. Exiting!\n", section);
                exit(EXIT_FAILURE);
            }
            else 
            {
                string_hash_insert(tmp->config_data.config, name, value);
            }
        }
    }
    else 
    {
        fprintf(stderr, "Invalid name %s in section %s on the config file!\n", name, section);
        return 0;
    }

    return 1;
}

void print_user_options (struct user_options *user_args)
{
    std::cout << "[main]" << std::endl;
    std::cout << "num_threads = " << user_args->num_threads << std::endl;
    std::cout << "dt = " << user_args->dt << std::endl;
    std::cout << "Final time = " << user_args->final_time << std::endl;
    std::cout << "Use Steady-State = " << user_args->use_steady_state << std::endl;
    std::cout << "Print rate = " << user_args->print_rate << std::endl;
    std::cout << "Steady-State rate = " << user_args->sst_rate << std::endl;
    std::cout << "Network filename = " << user_args->network_filename << std::endl;
    std::cout << "Steady-State filename = " << user_args->sst_filename << std::endl;
    std::cout << "Plot filename = " << user_args->plot_filename << std::endl << std::endl;
    
    std::cout << "[cell]" << std::endl;
    std::cout << "start_h = " << user_args->start_h << std::endl;
    std::cout << "num_div_cell = " << user_args->num_div_cell << std::endl;
    std::cout << "start diameter = " << user_args->start_diameter << std::endl;
    std::cout << "sigma_c = " << user_args->sigma_c << std::endl;
    std::cout << "G_gap = " << user_args->G_gap << std::endl;
    std::cout << "Model file path = " << user_args->model_file_path << std::endl << std::endl;

    if (user_args->stim_configs) 
    {
        std::cout << "[stimulus]" << std::endl;

        if (user_args->stim_configs->size == 1)
            printf ("Stimulus configuration:\n");
        else 
            printf ("Stimuli configuration:\n");

        for (int i = 0; i < user_args->stim_configs->size; i++) 
        {
            for (struct stim_config_elt *e = user_args->stim_configs->table[i % user_args->stim_configs->size]; e != 0;
                 e = e->next) {

                printf ("Stimulus name: %s\n", e->key);
                printf ("Stimulus start: %lf\n", e->value->stim_start);
                printf ("Stimulus duration: %lf\n", e->value->stim_duration);
                printf ("Stimulus current: %lf\n", e->value->stim_current);
                printf ("Start period: %lf\n", e->value->start_period);
                printf ("End period: %lf\n", e->value->end_period);
                printf ("Period step: %lf\n", e->value->period_step);
                printf ("Number of cycles: %d\n", e->value->n_cycles);
                printf ("Stimulus function: %s\n", e->value->config_data.function_name);
                struct string_hash *tmp = e->value->config_data.config;
                if (tmp->n == 1) 
                {
                    printf ("Stimulus extra parameter:\n");
                } 
                else if (tmp->n > 1) 
                {
                    printf ("Stimulus extra parameters:\n");
                }

                STRING_HASH_PRINT_KEY_VALUE_LOG (tmp);
            }
        }
    }

}

void display_usage (char *argv[])
{
    printf("=====================================================================================\n");
    printf("Usage:> %s -c <input_config_filename>\n",argv[0]);
    printf("-------------------------------------------------------------------------------------\n");
    printf("<input_config_filename> = Input filename with the parameters of the simulation\n");
    printf("=====================================================================================\n");
}