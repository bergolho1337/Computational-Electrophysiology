#ifndef CELULARMODEL_H
#define CELULARMODEL_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include "../include/model.h"
#include "../include/mitchell.h"
#include "../include/noble.h"
#include "../include/difrancesco.h"
#include "../include/luorudy.h"
#include "../include/timer.h"

using namespace std;

class CelullarModel
{
private:
  int ID;
  double DT;
  double TMAX;

  Model *MODEL;
  int NRATES;
  int NALGEBRAICS;
  int NCONSTANTS;
  double *STATES;
  double *RATES;
  double *CONSTANTS;
  double *ALGEBRAIC;
  
public:
  CelullarModel (int argc, char *argv[]);
  ~CelullarModel ();
  void buildMitchell ();
  void buildNoble ();
  void buildDiFrancesco ();
  void buildLuoRudy ();
  void Solve ();
  void PrintSolution (ofstream &out, const double t, const double STATES[]);
  void WriteLogFile (const double elapsedTime);
};

void Usage (const char pName[]);

#endif