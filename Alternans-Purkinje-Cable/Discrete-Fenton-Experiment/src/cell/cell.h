#ifndef _CELL_H_
#define _CELL_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>

struct cell_data
{
    double *yOld;
    double *yStar;
    double *yNew;
};

#endif