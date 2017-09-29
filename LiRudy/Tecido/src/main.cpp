/*
    Author: Lucas Berg
    ** Montar uma estrutura parecida com a do Noble para possibilitar a paralelizacao
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "../include/lirudy.h"
#include "../include/timer.h"

using namespace std;

/* Number of cells */
const int NC = 10;
/* Number of stimulus cells */
const int NSC = 2;

int main (int argc, char *argv[])
{
    if (argc-1 < 2)
    {
        cout << "-------- LiRudy 2011 -------------" << endl;
        cout << argv[0] << " <dt> <tmax>" << endl;
        cout << "----------------------------------" << endl;
        exit (EXIT_FAILURE);
    }
    LiRudy *lirudy = new LiRudy(argc,argv);
    //cout << *lirudy << endl;
    /*
    Cell *cells = (Cell*)malloc(sizeof(Cell)*NC);
    
    compGeometrics();

    setInitialConditions(cells,NC);
    //printCells(cells,NC);
    setStimulusCells(cells,NSC,NC);
    setTimeSettings(cells,NC);
    
    double start, finish, elapsed;
    GET_TIME(start);
    solveModel(cells,NC);
    GET_TIME(finish);
    elapsed = finish - start;
    printf("Time elapsed = %lf\n",elapsed);

    free(cells);
    */
    return 0;
}