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
#include "../include/fitzhugh.h"
#include "../include/tentusscherMC.h"
#include "../include/tentusscherEPI.h"
#include "../include/tentusscherENDO.h"

#include "../include/timer.h"

using namespace std;

class CelullarModel
{
private:
  int ID;
  int SOLVER;
  double DT;
  double TMAX;
  int PRINTRATE = 1;

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
  void buildFitzHugh ();
  void buildTenTusscherMC ();
  void buildTenTusscherEPI ();
  void buildTenTusscherENDO ();
  void Solve ();
  void ForwardEuler ();
  void RushLarsen ();
  
  void PrintSolution (ofstream &out, const int k, const double STATES[]);
  void WriteLogFile (const double elapsedTime);
};

void Usage (const char pName[]);

#endif