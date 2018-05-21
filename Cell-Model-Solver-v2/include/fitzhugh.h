// Impulses and physiological states in theoretical models of nerve membrane (1961)
// FitzHugh, R.A., 1961.

#ifndef FITZHUGH_H
#define FITZHUGH_H

#include "../include/model.h"

class FitzHugh: public Model
{
public:
    FitzHugh (string name) : Model(name) { };
    void initConst (double CONSTANTS[], double RATES[], double STATES[]);
    void compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    // Rush-Larsen
    void compRates_FE (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void RL (double DT, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
};


#endif