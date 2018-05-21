// A Model of Cardiac Electrical Activity Incorporating Ionic Pumps and Concentration Changes (1985)
// DiFrancesco, Noble, 1985

#ifndef DIFRANCESCO_H
#define DIFRANCESCO_H

#include "../include/model.h"

class DiFrancesco: public Model
{
public:
    DiFrancesco (string name) : Model(name) { };
    void initConst (double CONSTANTS[], double RATES[], double STATES[]);
    void compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    // Rush-Larsen
    void compRates_FE (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void RL (double DT, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
};


#endif