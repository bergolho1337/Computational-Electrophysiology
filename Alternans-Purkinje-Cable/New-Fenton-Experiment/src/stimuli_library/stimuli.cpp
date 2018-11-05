//
// Created by sachetto on 13/10/17.
//

#include <cstdio>
#include <cstdbool>
#include <cstdint>
#include <cstddef>


//#include "../utils/utils.h"
#include "../monodomain/constants.h"
//#include "../alg/grid/grid.h"
#include "../monodomain/config/stim_config.h"
//#include "../libraries_common/config_helpers.h"

extern "C" SET_SPATIAL_STIM(stim_if_id_less_than) 
{
/*
    uint32_t n_active = the_grid->num_active_cells;
    //struct cell_node **ac = the_grid->active_cells;

    bool stim;
    real stim_current = config->stim_current;
    real stim_value;

    int id;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, id, config->config_data.config, "id_limit");

    if(config->spatial_stim_currents) 
    {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

	int i;

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) 
    {
        stim = i <= id;

        if (stim) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
    */
}

extern "C" SET_SPATIAL_STIM(stim_if_id_greater_than) {

    /*
    uint32_t n_active = the_grid->num_active_cells;
    //struct cell_node **ac = the_grid->active_cells;

    bool stim;
    real stim_current = config->stim_current;
    real stim_value;

    int id;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, id, config->config_data.config, "id_limit");

    if(config->spatial_stim_currents) 
    {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

	int i;

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) 
    {
        stim = i >= id;

        if (stim) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
    */
}