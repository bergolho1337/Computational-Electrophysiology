// Alternans and spiral breakup in a human ventricular tissue model
// Endocardium cell
// TenTusscher, 2006

#ifndef TENTUSSCHERENDO_H
#define TENTUSSCHERENDO_H

#include "../include/model.h"

class TenTusscherENDO: public Model
{
public:
    TenTusscherENDO (string name) : Model(name) { };
    void initConst (double CONSTANTS[], double RATES[], double STATES[]);
    void compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    // Rush-Larsen
    void compRates_FE (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void RL (double DT, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
};

#endif