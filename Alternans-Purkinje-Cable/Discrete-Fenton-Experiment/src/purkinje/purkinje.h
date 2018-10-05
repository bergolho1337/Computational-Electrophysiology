#ifndef _PURKINJE_H_
#define _PURKINJE_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "../graph/graph.h"
#include "../config/user_config.h"

using namespace std;

void configure_purkinje_network_from_options (struct graph *the_purkinje_network, struct user_options *config);
void set_purkinje_network_from_file(struct graph *the_purkinje_network, char *network_filename, const double size_volume, const int num_div_cell);
void build_skeleton_mesh (char *filename, struct graph *sk);
void build_mesh_purkinje (struct graph *the_purkinje_network, struct graph *the_skeleton, const double size_volume, const int num_div_cell);

#endif