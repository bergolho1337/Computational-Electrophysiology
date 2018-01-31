#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;

// Abstract class
class Model
{
protected:
  string NAME;
public:
  Model (string name) : NAME(name) { };
  // Abstract functions
  virtual void initConst (double CONSTANTS[], double RATES[], double STATES[]) = 0;
  virtual void compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]) = 0;
  virtual void compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[]) = 0;
};

#endif