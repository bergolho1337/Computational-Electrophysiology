#include "grid.h"

struct grid* new_grid()
{
    struct grid* result = (struct grid*) malloc(sizeof(struct grid));
    result->first_cell = NULL;
    result->active_cells = NULL;

    // Purkinje
    result->the_purkinje_network = new_graph();

    return result;
}