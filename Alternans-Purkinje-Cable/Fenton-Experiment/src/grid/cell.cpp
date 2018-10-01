#include "cell.h"

struct cell_node* new_cell_node() 
{
    struct cell_node* result = (struct cell_node*)malloc(sizeof(struct cell_node));
    init_cell_node(result);
    return result;

}

void init_cell_node(struct cell_node *cell_node) 
{

    cell_node->center_x = 0.0;
    cell_node->center_y = 0.0;
    cell_node->center_z = 0.0;

    cell_node->active = true;

    cell_node->previous    = NULL;
    cell_node->next        = NULL;

    cell_node->grid_position        = 0;
    cell_node->sv_position          = 0;
    cell_node->face_length          = 1.0;
    cell_node->half_face_length     = 0.5;

    cell_node->v = 0;

    cell_node->b = 0.0;

}

void set_cell_node_data(struct cell_node *the_cell, double face_length, double half_face_length,
                                void *previous, void *next,
                                uint32_t grid_position,
                                double center_x, double center_y, double center_z)
{
    the_cell->face_length = face_length;
    the_cell->half_face_length = half_face_length;
    the_cell->previous = (struct cell_node*)previous;
    the_cell->next = (struct cell_node*)next;
    the_cell->grid_position = grid_position;
    the_cell->center_x = center_x;
    the_cell->center_y = center_y;
    the_cell->center_z = center_z;
}

void free_cell_node(struct cell_node *cell_node) 
{
    free(cell_node);
}

