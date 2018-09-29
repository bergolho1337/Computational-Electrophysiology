//
// Created by sachetto on 13/10/17.
// Edited by bergolho on 16/05/18
//

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "../utils/utils.h"
#include "../monodomain/constants.h"
#include "../alg/grid/grid.h"
#include "../monodomain/config/assembly_matrix_config.h"
#include "../libraries_common/config_helpers.h"

static inline double ALPHA (double beta, double cm, double dt, double h) {
    return (((beta * cm) / dt) * UM2_TO_CM2) * pow (h, 3.0);
}

void initialize_diagonal_elements (struct monodomain_solver *the_solver, struct grid *the_grid) {

    double alpha, h;
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    double beta = the_solver->beta;
    double cm = the_solver->cm;

    double dt = the_solver->dt;

    int i;

#pragma omp parallel for private(alpha, h)
    for (i = 0; i < num_active_cells; i++) {
        h = ac[i]->face_length;
        alpha = ALPHA (beta, cm, dt, h);

        struct element element;
        element.column = ac[i]->grid_position;
        element.cell = ac[i];
        element.value = alpha;

        if (ac[i]->elements != NULL) {
                    sb_free (ac[i]->elements);
        }

        ac[i]->elements = NULL;
                sb_reserve (ac[i]->elements, 7);
                sb_push (ac[i]->elements, element);
    }
}


static void fill_discretization_matrix_elements (double sigma_x, double sigma_y, double sigma_z, struct cell_node *grid_cell,
                                                 void *neighbour_grid_cell, char direction) {

    uint32_t position;
    bool has_found;
    double h;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    double sigma_x1 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);
    double sigma_x2 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);
    double sigma_y1 = (2.0f * sigma_y * sigma_y) / (sigma_y + sigma_y);
    double sigma_y2 = (2.0f * sigma_y * sigma_y) / (sigma_y + sigma_y);
    double sigma_z1 = (2.0f * sigma_z * sigma_z) / (sigma_z + sigma_z);
    double sigma_z2 = (2.0f * sigma_z * sigma_z) / (sigma_z + sigma_z);

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
     * cell which is a cell node. */
    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    char neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    if (neighbour_grid_cell_level > grid_cell->cell_data.level) {
        if ((neighbour_grid_cell_type == TRANSITION_NODE_TYPE)) {
            has_found = false;
            while (!has_found) {
                if (neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {
                    white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
                    if (white_neighbor_cell->single_connector == NULL) {
                        has_found = true;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } else {
                    break;
                }
            }
        }
    } else {
        if (neighbour_grid_cell_level <= grid_cell->cell_data.level &&
            (neighbour_grid_cell_type == TRANSITION_NODE_TYPE)) {
            has_found = false;
            while (!has_found) {
                if (neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {
                    white_neighbor_cell = (struct transition_node *)(neighbour_grid_cell);
                    if (white_neighbor_cell->single_connector == 0) {
                        has_found = true;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } else {
                    break;
                }
            }
        }
    }

    // Tratamos somente os pontos interiores da malha.
    if (neighbour_grid_cell_type == CELL_NODE_TYPE) {

        black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);

        if (black_neighbor_cell->active) {

            if (black_neighbor_cell->cell_data.level > grid_cell->cell_data.level) {
                h = black_neighbor_cell->face_length;
            } else {
                h = grid_cell->face_length;
            }

            lock_cell_node (grid_cell);

            struct element *cell_elements = grid_cell->elements;
            position = black_neighbor_cell->grid_position;

            size_t max_elements = sb_count (cell_elements);
            bool insert = true;

            for (size_t i = 1; i < max_elements; i++) {
                if (cell_elements[i].column == position) {
                    insert = false;
                    break;
                }
            }

            // TODO: maybe each element can have a different sigma!
            if (insert) {

                struct element new_element;
                new_element.column = position;
                if (direction == 'n') { // Z direction
                    new_element.value = -sigma_z1 * h;
                    cell_elements[0].value += (sigma_z1 * h);
                } else if (direction == 's') { // Z direction
                    new_element.value = -sigma_z2 * h;
                    cell_elements[0].value += sigma_z2 * h;
                } else if (direction == 'e') { // Y direction
                    new_element.value = -sigma_y1 * h;
                    cell_elements[0].value += (sigma_y1 * h);
                } else if (direction == 'w') { // Y direction
                    new_element.value = -sigma_y2 * h;
                    cell_elements[0].value += (sigma_y2 * h);
                } else if (direction == 'f') { // X direction
                    new_element.value = -sigma_x1 * h;
                    cell_elements[0].value += (sigma_x1 * h);
                } else if (direction == 'b') { // X direction
                    new_element.value = -sigma_x2 * h;
                    cell_elements[0].value += (sigma_x2 * h);
                }

                new_element.cell = black_neighbor_cell;
                        sb_push (grid_cell->elements, new_element);
            }
            unlock_cell_node (grid_cell);

            lock_cell_node (black_neighbor_cell);
            cell_elements = black_neighbor_cell->elements;
            position = grid_cell->grid_position;

            max_elements = sb_count (cell_elements);

            insert = true;
            for (size_t i = 1; i < max_elements; i++) {
                if (cell_elements[i].column == position) {
                    insert = false;
                    break;
                }
            }

            if (insert) {

                struct element new_element;
                new_element.column = position;
                if (direction == 'n') { // Z direction
                    new_element.value = -sigma_z1 * h;
                    cell_elements[0].value += (sigma_z1 * h);
                } else if (direction == 's') { // Z direction
                    new_element.value = -sigma_z2 * h;
                    cell_elements[0].value += (sigma_z2 * h);
                } else if (direction == 'e') { // Y direction
                    new_element.value = -sigma_y1 * h;
                    cell_elements[0].value += (sigma_y1 * h);
                } else if (direction == 'w') { // Y direction
                    new_element.value = -sigma_y2 * h;
                    cell_elements[0].value += (sigma_y2 * h);
                } else if (direction == 'f') { // X direction
                    new_element.value = -sigma_x1 * h;
                    cell_elements[0].value += (sigma_x1 * h);
                } else if (direction == 'b') { // X direction
                    new_element.value = -sigma_x2 * h;
                    cell_elements[0].value += (sigma_x2 * h);
                }

                new_element.cell = grid_cell;
                        sb_push (black_neighbor_cell->elements, new_element);
            }

            unlock_cell_node (black_neighbor_cell);
        }
    }
}

void initialize_diagonal_elements_purkinje (struct monodomain_solver *the_solver, struct grid *the_grid) 
{

    double alpha, h;
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    struct node *n = the_grid->the_purkinje_network->list_nodes;
    double beta = the_solver->beta;
    double cm = the_solver->cm;

    double dt = the_solver->dt;

    int i;

    for (i = 0; i < num_active_cells; i++) 
    {
        h = ac[i]->face_length;
        alpha = ALPHA (beta, cm, dt, h);

        struct element element;
        element.column = ac[i]->grid_position;
        element.cell = ac[i];
        element.value = alpha;

        if (ac[i]->elements != NULL) 
        {
            sb_free (ac[i]->elements);
        }

        ac[i]->elements = NULL;
        sb_reserve (ac[i]->elements,n->num_edges);
        sb_push (ac[i]->elements, element);

        n = n->next;
    }       
}

static void fill_discretization_matrix_elements_purkinje (double sigma_x, struct cell_node **grid_cells, uint32_t num_active_cells,
                                                        struct node *pk_node) 
{

    //uint32_t position;
    //bool has_found;
    double h;

    double sigma_x1 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);
    //double sigma_x2 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);

    struct edge *e;
    struct element **cell_elements;

    int i;

    for (i = 0; i < num_active_cells; i++, pk_node = pk_node->next)
    {
        cell_elements = &grid_cells[i]->elements;
        h = grid_cells[i]->face_length;

        e = pk_node->list_edges;

        //printf("\n");
        // Do the mapping of the edges from the graph to the sparse matrix data structure ...
        while (e != NULL)
        {
            struct element new_element;

            // Neighbour elements ...
            new_element.column = e->id;
            new_element.value = -sigma_x1 * h;
            new_element.cell = grid_cells[e->id];

            // Diagonal element ...
            cell_elements[0]->value += (sigma_x1 * h);
            //printf("Node %d -- Diagonal = %lf\n",pk_node->id,cell_elements[0]->value);

            sb_push(grid_cells[i]->elements,new_element);   

            e = e->next;         
        }
    }
}

ASSEMBLY_MATRIX(no_fibers_assembly_matrix) 
{

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements (the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");

    #pragma omp parallel for
        for (i = 0; i < num_active_cells; i++) 
        {

            // Computes and designates the flux due to south cells.
            fill_discretization_matrix_elements (sigma_x, sigma_y, sigma_z, ac[i], ac[i]->south, 's');

            // Computes and designates the flux due to north cells.
            fill_discretization_matrix_elements (sigma_x, sigma_y, sigma_z, ac[i], ac[i]->north, 'n');

            // Computes and designates the flux due to east cells.
            fill_discretization_matrix_elements (sigma_x, sigma_y, sigma_z, ac[i], ac[i]->east, 'e');

            // Computes and designates the flux due to west cells.
            fill_discretization_matrix_elements (sigma_x, sigma_y, sigma_z, ac[i], ac[i]->west, 'w');

            // Computes and designates the flux due to front cells.
            fill_discretization_matrix_elements (sigma_x, sigma_y, sigma_z, ac[i], ac[i]->front, 'f');

            // Computes and designates the flux due to back cells.
            fill_discretization_matrix_elements (sigma_x, sigma_y, sigma_z, ac[i], ac[i]->back, 'b');
        }

    /*
    for (i = 0; i < num_active_cells; i++)
    {
        printf("\nCell %d -- Diagonal = %lf\n",i,ac[i]->elements[0].value);
        int count = sb_count(ac[i]->elements);
        printf("\tElements:\n");
        for (int j = 1; j < count; j++)
            printf("\t%d -- Column = %d -- Value = %lf\n",ac[i]->elements[j].column,ac[i]->elements[j].column,ac[i]->elements[j].value);
    }

    printf("Leaving program ...\n");
    exit(EXIT_FAILURE);
    */
}

ASSEMBLY_MATRIX(purkinje_fibers_assembly_matrix) 
{

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    struct node *pk_node = the_grid->the_purkinje_network->list_nodes;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    fill_discretization_matrix_elements_purkinje(sigma_x,ac,num_active_cells,pk_node);

    /*
    for (int i = 0; i < num_active_cells; i++)
    {
        printf("\nCell %d -- Diagonal = %lf\n",i,ac[i]->elements[0].value);
        int count = sb_count(ac[i]->elements);
        printf("\tElements:\n");
        for (int j = 1; j < count; j++)
            printf("\t%d -- Column = %d -- Value = %lf\n",ac[i]->elements[j].column,ac[i]->elements[j].column,ac[i]->elements[j].value);
    }

    printf("Leaving program ...\n");
    exit(EXIT_FAILURE);
    */

}