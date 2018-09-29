#include "purkinje_config.h"

struct purkinje_config* new_purkinje_config ()
{
    struct purkinje_config *config = (struct purkinje_config*)malloc(sizeof(struct purkinje_config));

    config->network_filename = NULL;
    config->network_filename_was_set = false;

    config->name = NULL;
    config->name_was_set = false;

    config->start_h = 100.0;
    config->start_h_was_set = false;

    return config;
}