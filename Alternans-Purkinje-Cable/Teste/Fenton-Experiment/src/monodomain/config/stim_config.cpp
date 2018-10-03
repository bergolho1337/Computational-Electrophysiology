//
// Created by sachetto on 13/10/17.
//

#include "stim_config.h"
#include <dlfcn.h>
#include <cstring>

void init_stim_functions(struct stim_config *config, char* stim_name) {

}

struct stim_config* new_stim_config() {
    struct stim_config *result = (struct stim_config*) malloc(sizeof(struct stim_config));
    init_config_common_data(&(result->config_data));
    result->spatial_stim_currents = NULL;

    result->stim_current_was_set = false;
    result->stim_duration_was_set = false;
    result->stim_start_was_set = false;
    return result;
}

void print_stim_config_values(struct stim_config* s) {
    fprintf(stdout,"stim_start %.10lf\n",s->stim_start);
    fprintf(stdout,"stim_duration %.10lf\n",s->stim_duration);
    fprintf(stdout,"stim_current %.10lf\n",s->stim_current);
    fprintf(stdout,"start_period %.10lf\n",s->start_period);
    fprintf(stdout,"end_period %.10lf\n",s->end_period);
    fprintf(stdout,"period_step %.10lf\n",s->period_step);
    fprintf(stdout,"n_cycles %d\n",s->n_cycles);
    fprintf(stdout,"id_limit %d\n",s->id_limit);
    //print_to_stdout_and_file("stim_duration %lf\n", s->stim_duration);
    //print_to_stdout_and_file("stim_current %lf\n", s->stim_current);
    //print_to_stdout_and_file("stim_function %s\n", s->config_data.function_name);
}

void free_stim_config(struct stim_config* s) {
    free(s->config_data.function_name);
    free(s->config_data.library_file_path);
    free(s->spatial_stim_currents);
    if(s->config_data.handle)
        dlclose(s->config_data.handle);
    //string_hash_destroy(s->config_data.config);
    free(s);
}
