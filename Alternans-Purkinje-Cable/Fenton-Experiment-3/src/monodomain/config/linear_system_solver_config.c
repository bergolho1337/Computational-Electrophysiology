//
// Created by sachetto on 13/10/17.
// Edited by bergolho on 16/05/18
//

#include "linear_system_solver_config.h"

#ifdef _MSC_VER
#include "../../dlfcn-win32/dlfcn.h"
#else
#include <dlfcn.h>
#endif
#include <string.h>
#include "../../utils/logfile_utils.h"

void init_linear_system_solver_functions(struct linear_system_solver_config *config) 
{

    char *function_name = config->config_data.function_name;

#ifdef _MSC_VER
	char *default_function = "./shared_libs/libdefault_linear_system_solver.dll";
#else
	char *default_function = "./shared_libs/libdefault_linear_system_solver.so";
#endif

    if(config->config_data.library_file_path == NULL) 
    {
        print_to_stdout_and_file("Using the default library to solve the linear system\n");
        config->config_data.library_file_path = strdup(default_function);
        config->config_data.library_file_path_was_set = true;
    }
    else 
    {
        print_to_stdout_and_file("Opening %s as linear system solver lib\n", config->config_data.library_file_path);
    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) 
    {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name)
    {
        config->solve_linear_system = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL)  
        {
            fprintf(stderr, "\n%s function not found in the provided linear_system_solver library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else 
    {
        fprintf(stderr, "No function name for matrix assembly library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct linear_system_solver_config* new_linear_system_solver_config() 
{
    struct linear_system_solver_config *result = (struct linear_system_solver_config*) malloc(sizeof(struct linear_system_solver_config));

    init_config_common_data(&(result->config_data));
    result->solve_linear_system = NULL;
    return result;
}

void print_linear_system_solver_config_values(struct linear_system_solver_config* s) 
{
    printf("linear_system_solver_function: %s\n",s->config_data.function_name);
    printf("linear_system_solver_library_file: %s\n",s->config_data.library_file_path);
    printf("linear_system_solver_config:\n");
    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_linear_system_solver_config(struct linear_system_solver_config* s) 
{

    free(s->config_data.library_file_path);
    free(s->config_data.function_name);
    string_hash_destroy(s->config_data.config);

    if(s->config_data.handle)
        dlclose(s->config_data.handle);
    free(s);
}
