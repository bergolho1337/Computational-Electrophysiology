#ifndef _CELL_H_
#define _CELL_H_

#include <cstdlib>
#include <cstdio>
#include <cstdbool>
#include <cstdint>

struct element 
{
    double value;
    uint32_t column; // Column of the matrix to which this element belongs.
    struct cell_node *cell;
};

struct cell_node 
{
    bool active;

    // Cell position
    double center_x, center_y, center_z;

    struct cell_node *previous; // Previous cell.
    struct cell_node *next;     // Next cell.

    // Indicates position of cell on grid
    uint32_t grid_position;

    // Cell geometry.
    double half_face_length;
    double diameter;

    //______________________________________________________________________________
    /* The matrix row. The elements[0] corresponds to the diagonal element of the row. */
    //element_vector *elements;
    //struct element *elements;

    //______________________________________________________________________________
    /* Variables used in solving the discretized system Ax = b through the 
    LU Decomposition */
    double b;

    // Variables used by some applications of partial differential equations.
    double v;

    uint32_t sv_position;
    double face_length;

};

struct cell_node* new_cell_node();
void free_cell_node(struct cell_node *cell_node);
void init_cell_node(struct cell_node *cell_node);
void set_cell_node_data(struct cell_node *the_cell, double face_length, double half_face_length,
                                void *previous, void *next,
                                uint32_t grid_position,
                                double center_x, double center_y, double center_z); 


#endif