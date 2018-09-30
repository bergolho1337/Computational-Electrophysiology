//
// Created by bergolho on 19/07/18.
//

#include "purkinje_helpers.h"

#include "../libraries_common/config_helpers.h"
#include "../monodomain/config/purkinje_config.h"
#include "../utils/logfile_utils.h"
#include <assert.h>
#include <time.h>

#ifdef _MSC_VER
    #include <process.h>
    #define getpid _getpid
#else
    #include <unistd.h>
#endif

SET_SPATIAL_PURKINJE (initialize_purkinje_with_custom_mesh) 
{

    char *network_file;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR (network_file,config->config_data.config,"network_file");

    // TO DO: Consider also the diameter of the Purkinje cell ...
    double side_length = config->start_h;
    
    print_to_stdout_and_file("Loading a custom Purkinje Network: %s\n",config->domain_name);
    print_to_stdout_and_file("Using the Purkinje library function: \"initialize_purkinje_with_custom_mesh\"\n");
    print_to_stdout_and_file("Discretization for the Purkinje Network Mesh: %.10lf um\n",side_length);
    set_custom_purkinje_network(the_grid, network_file, side_length);

    // TO DO: Populate the 'grid_cell' linked-list with the nodes from the graph
    //        Some parameters from the 'cell_node' structure will not be used
    initialize_and_construct_grid_purkinje(the_grid);

    free (network_file);

}

// TO DO: Build some benchmark Purkinje network for quick tests ...