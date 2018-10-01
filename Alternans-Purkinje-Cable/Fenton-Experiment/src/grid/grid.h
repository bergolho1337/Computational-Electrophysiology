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

void initialize_and_construct_grid_purkinje (struct grid *grid);
void initialize_grid_purkinje (struct grid *grid);
void construct_grid_purkinje (struct grid *the_grid);

void order_grid_cells (struct grid *the_grid);
void save_old_cell_positions (struct grid *the_grid);

bool print_grid_and_check_for_activity(const struct grid *the_grid, FILE *output_file, const int count, const bool binary);

#endif