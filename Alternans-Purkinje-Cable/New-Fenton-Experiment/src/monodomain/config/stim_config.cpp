//
// Created by sachetto on 13/10/17.
//

#include "stim_config.h"
#include <dlfcn.h>
#include <cstring>
//#include "../../utils/logfile_utils.h"
/*
void init_stim_functions(struct stim_config *config, char* stim_name) {

#ifdef _MSC_VER
	char *default_library_name = "./shared_libs/default_stimuli.dll";
#else
	char *default_library_name = "./shared_libs/libdefault_stimuli.so";
#endif

    char *function_name = config->config_data.function_name;

    if(config->config_data.library_file_path == NULL) {
        printf("Using the default library for stimuli functions for %s\n", stim_name);
        config->config_data.library_file_path = strdup(default_library_name);
        config->config_data.library_file_path_was_set = true;
    }
    else {
        printf("Opening %s as stimuli lib for %s\n", config->config_data.library_file_path, stim_name);
    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name) {
        config->set_spatial_stim = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL) {
            fprintf(stderr, "\n%s function not found in the provided stimuli library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for sitimuli library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }


}
*/

struct stim_config* new_stim_config() {
    struct stim_config *result = (struct stim_config*) malloc(sizeof(struct stim_config));
    init_config_common_data(&(result->config_data));
    //result->set_spatial_stim = NULL;
    result->spatial_stim_currents = NULL;

    result->stim_current_was_set = false;
    result->stim_duration_was_set = false;
    result->stim_start_was_set = false;
    return result;
}

void print_stim_config_values(struct stim_config* s) {
    printf("stim_start %lf\n", s->stim_start);
    printf("stim_duration %lf\n", s->stim_duration);
    printf("stim_current %lf\n", s->stim_current);
    printf("start_period %lf\n", s->start_period);
    printf("end_period %lf\n", s->end_period);
    printf("period_step %lf\n", s->period_step);
    printf("n_cycles %d\n", s->n_cycles);
    //printf("stim_function %s\n", s->config_data.function_name);
}

void free_stim_config(struct stim_config* s) {
    free(s->config_data.function_name);
    free(s->config_data.library_file_path);
    free(s->spatial_stim_currents);
    if(s->config_data.handle)
        dlclose(s->config_data.handle);
    string_hash_destroy(s->config_data.config);
    free(s);
}