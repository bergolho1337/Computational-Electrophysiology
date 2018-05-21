// Alternans and spiral breakup in a human ventricular tissue model
// Midmiocardium cell
// TenTusscher, 2006

#ifndef TENTUSSCHERMC_H
#define TENTUSSCHERMC_H

#include "../include/model.h"

class TenTusscherMC: public Model
{
public:
    TenTusscherMC (string name) : Model(name) { };
    void initConst (double CONSTANTS[], double RATES[], double STATES[]);
    void compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    // Rush-Larsen
    void compRates_FE (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void RL (double DT, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
};

#endif