//
// Created by bergolho on 19/07/18.
//

#include "purkinje_config.h"
#include <dlfcn.h>
#include <cstring>
#include "../../utils/logfile_utils.h"

void init_purkinje_functions (struct purkinje_config *config) 
{

}

struct purkinje_config* new_purkinje_config() 
{
    struct purkinje_config *result = (struct purkinje_config*) malloc(sizeof(struct purkinje_config));

    init_config_common_data(&(result->config_data));

    result->start_h_was_set = false;
    result->network_name_was_set = false;

    result->set_spatial_purkinje = NULL;
    result->network_name = NULL;
    result->network_filename = NULL;
    return result;
}

void print_purkinje_config_values (struct purkinje_config* s) 
{
    printf("purkinje_name: %s\n",s->network_name);
    printf("purkinje_filename: %s\n",s->network_filename);
    printf("start_discretization: %lf\n",s->start_h);
}

void free_purkinje_config(struct purkinje_config* s) 
{
    //free(s->config_data.library_file_path);
    //free(s->config_data.function_name);
    //string_hash_destroy(s->config_data.config);
    free(s->network_name);
    free(s->network_filename);
    if(s->config_data.handle)
        dlclose(s->config_data.handle);
    free(s);
}
