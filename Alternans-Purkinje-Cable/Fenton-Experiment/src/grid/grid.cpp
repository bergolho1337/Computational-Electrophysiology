#include "grid.h"

struct grid* new_grid()
{
    struct grid *result = (struct grid*)malloc(sizeof(struct grid));

    result->the_purkinje_network = NULL;
    
    return result;
}

void configure_grid_from_options (struct grid *the_grid, struct user_options *options)
{
    char* network_filename = options->network_filename;

    the_grid->dx = options->start_h / options->num_div_cell;

    the_grid->the_purkinje_network = set_purkinje_mesh_from_file(network_filename,the_grid->dx);
    //the_grid->the_purkinje_network->set_gap_junctions(options->num_div_cell);
    the_grid->the_purkinje_network->set_no_gap_junctions();
    //the_grid->the_purkinje_network->print();

}

void clean_and_free_grid(struct grid* the_grid)
{
    assert(the_grid);

    clean_grid(the_grid);

    free (the_grid);
}

void clean_grid(struct grid *the_grid)
{
    assert(the_grid);
    
    if (the_grid->the_purkinje_network) 
    {
        delete the_grid->the_purkinje_network;
    }

}