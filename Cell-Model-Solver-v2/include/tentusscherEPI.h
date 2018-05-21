// Alternans and spiral breakup in a human ventricular tissue model
// Epicardium cell
// TenTusscher, 2006

#ifndef TENTUSSCHEREPI_H
#define TENTUSSCHEREPI_H

#include "../include/model.h"

class TenTusscherEPI: public Model
{
public:
    TenTusscherEPI (string name) : Model(name) { };
    void initConst (double CONSTANTS[], double RATES[], double STATES[]);
    void compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    // Rush-Larsen
    void compRates_FE (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void RL (double DT, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
};

#endif