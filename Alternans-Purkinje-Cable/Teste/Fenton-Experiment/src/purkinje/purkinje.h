#ifndef _PURKINJE_H_
#define _PURKINJE_H_

#include "../monodomain/config/purkinje_config.h"
#include "../grid/grid.h"

void set_spatial_purkinje (struct purkinje_config *pk_config, struct grid *grid);
void set_custom_purkinje_network (struct grid *the_grid, const char *file_name,\
                                    const double side_length, const double diameter);
void set_purkinje_network_from_file (struct graph *the_purkinje_network, const char *file_name,\
                                        const double side_length);
void build_skeleton_purkinje (const char *filename, struct graph *skeleton_network);
void build_mesh_purkinje (struct graph *the_purkinje_network, struct graph *skeleton_network, const double side_length);
void grow_segment (struct graph *the_purkinje_network, struct node *u, struct edge *v, uint32_t *map_skeleton_to_mesh);
void depth_first_search (struct graph *the_purkinje_network, struct node *u, int level, uint32_t *map_skeleton_to_mesh);
void calc_unitary_vector (double d_ori[], struct node *u, struct node *v);
void write_purkinje_network_to_vtk (struct graph *the_purkinje_network);

#endif