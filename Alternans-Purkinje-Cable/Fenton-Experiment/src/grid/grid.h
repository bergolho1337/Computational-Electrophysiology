#ifndef _GRID_H_
#define _GRID_H_

#include <cstdlib>
#include <cstdio>
#include "cell.h"
#include "../graph/graph.h"

struct grid 
{

    struct cell_node *first_cell;     // First cell of grid.
    double side_length;        // Length of cube grid. Default = 1.0.
    uint32_t number_of_cells;  // Number of cells of grid.

    uint32_t num_active_cells;

    struct cell_node* *active_cells;

    // Purkinje network graph
    struct graph *the_purkinje_network;
};

struct grid* new_grid();

#endif