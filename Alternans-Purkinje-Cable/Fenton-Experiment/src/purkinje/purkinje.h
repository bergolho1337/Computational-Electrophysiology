#ifndef FENTON_PURKINJE_H_
#define FENTON_PURKINJE_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <algorithm>

#include "graph.h"

using namespace std;

Graph* set_purkinje_mesh_from_file (const string network_filename, const double h);
Graph* initialize_purkinje_mesh (const string filename, const double h);
void build_skeleton_purkinje (Graph *sk, const string filename);
void build_mesh_purkinje (Graph *mesh, Graph *sk, const double h);
void grow_segment (Graph *mesh, Node *u, Edge *v, int *map_skeleton_to_mesh, const double h);
void depth_first_search(Graph *mesh, Node *u, int level, int *map_skeleton_to_mesh, const double h);
void calc_unitary_vector (double d_ori[], Node *u, Node *v);

#endif