// A Model of the Ventricular Cardiac Action Potential - Depolarisation, Repolarisation and Their Interaction, Ching-hsing Luo and Yoram Rudy, 1991
// Luo Rudy 1991

#ifndef LUORUDY_H
#define LUORUDY_H

#include "../include/model.h"

class LuoRudy: public Model
{
public:
    LuoRudy (string name) : Model(name) { };
    void initConst (double CONSTANTS[], double RATES[], double STATES[]);
    void compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    // Rush-Larsen
    void compRates_FE (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
    void RL (double DT, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]);
};


#endif