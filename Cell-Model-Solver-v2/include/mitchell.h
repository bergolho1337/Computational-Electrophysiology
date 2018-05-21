// A two-current model for the dynamics of cardiac membrane (2003)
// Mitchell, Schaeffer, 2003

#ifndef MITCHELL_H
#define MITCHELL_H

#include "../include/model.h"

class Mitchell: public Model
{
public:
    Mitchell (string name) : Model(name) { };
    void initConst (double CONSTANTS[], double RATES[], double STATES[]);
    void compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    // Rush-Larsen
    void compRates_FE (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void RL (double DT, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
};


#endif