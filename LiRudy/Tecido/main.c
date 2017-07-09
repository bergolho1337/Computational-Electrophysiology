/*
    Author: Lucas Berg
*/

#include <stdio.h>
#include <stdlib.h>
#include "prd.h"

/* Number of cells */
const int NC = 10;
/* Number of stimulus cells */
const int NSC = 2;

int main ()
{
    Cell *cells = (Cell*)malloc(sizeof(Cell)*NC);
    
    compGeometrics();

    setInitialConditions(cells,NC);
    //printCells(cells,NC);
    setStimulusCells(cells,NSC,NC);
    setTimeSettings(cells,NC);
    
    solveModel(cells,NC);

    free(cells);
    return 0;
}