// A Modification of the Hodgkin-Huxley Equations Applicable to Purkinje Fibre Action and Pace-Maker Potentials
// Noble, 1962

#ifndef NOBLE_H
#define NOBLE_H

#include "../include/model.h"

class Noble: public Model
{
public:
    Noble (string name) : Model(name) { };
    void initConst (double CONSTANTS[], double RATES[], double STATES[]);
    void compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
};


#endif