//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_GRID_H
#define MONOALG3D_GRID_H

#include "../cell/cell.h"
#include "../../graph/graph.h"

#include <stdlib.h>
#include <stdio.h>

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
void clean_and_free_grid(struct grid* the_grid);

void print_grid(struct grid* the_grid, FILE *output_file);

void clean_grid(struct grid *the_grid);
void order_grid_cells (struct grid *the_grid);

void set_grid_flux(struct grid *the_grid);

void print_grid_matrix(struct grid *the_grid, FILE* output_file);
void print_grid_vector(struct grid* the_grid, FILE *output_file, char name);
double * grid_vector_to_array(struct grid *the_grid, char name, uint32_t *num_lines);

void save_grid_domain (struct grid * the_grid, const char *file_name);

void lock_grid(struct grid *the_grid);

void unlock_grid(struct grid *the_grid);


void initialize_and_construct_grid_purkinje (struct grid *the_grid);
void initialize_grid_purkinje (struct grid *the_grid);
void construct_grid_purkinje (struct grid *the_grid);

#endif //MONOALG3D_GRID_H
