#ifndef NOBLE_H_
#define NOBLE_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <omp.h>
#include "../include/cell.h"

using namespace std;

class Noble
{
    static constexpr int NCELLS = 500;
    static constexpr int NEQ = 4;
    static constexpr int OUT_ID = 4;

private:
    int nthreads;
    int nbeats;
    int n;
    double dt;
    double tmax;
    vector<Cell> cells;
public:
    Noble (int argc, char *argv[]);
    void allocMem ();
    void setInitCond ();
    void solve ();

    // DEBUG
    friend ostream& operator<< (ostream &ost, const Noble &noble)
    {
        ost << "-------- PARAMETERS ------------" << endl;
        ost << "Number of beats = " << noble.nbeats << endl;
        ost << "Dt = " << noble.dt << endl;
        ost << "tmax = " << noble.tmax << endl;
        ost << "Number of timesteps = " << noble.n << endl;
        ost << "--------------------------------" << endl;
        for (int i = 0; i < (int)noble.cells.size(); i++)
        {
            ost << "Cell " << i << endl;
            ost << noble.cells[i] << endl;
        }
        ost << "--------------------------------" << endl;  

        return ost;
    }
};

void readInput (int argc, char *argv[]);
void setInitialConditions (double *yOld);
void writeIteration (double t, double *y);
void solveEDO ();
void isEquilibrium (double t, double *yOld, double *yNew);

#endif
