#include "monodomain_solver.h"

struct monodomain_solver *new_monodomain_solver ()
{
    struct monodomain_solver *result = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));

    result->beta = 0.14;
    result->cm = 1.2;

    return result;
}
