#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "mitchell.h"
#include "noble.h"
#include "hodkin.h"
#include "fitz.h"

// Function pointer
typedef double (*Func) (int k, double t, double *y);

struct CelularModel
{
  int id;                           // Id of the celular model choosen by the user
  int num_eq;                       // Number of equations
  double *y0;                       // Vector of the initial conditions
  double dt;                        // Timestep of the solver
  double t_max;                     // Max time of the simulation
  Func *f;                          // Vector of the functions of the method
}typedef CelularModel;

CelularModel* initModel (int argc, char *argv[]);

void buildMitchell (CelularModel *cm);
void buildNoble (CelularModel *cm);
void buildHodkin (CelularModel *cm);
void buildFitzhugh (CelularModel *cm);
void buildLiRudy (CelularModel *cm);
void buildHaq (CelularModel *cm);

void solveModel (CelularModel *cm);
void freeModel (CelularModel *cm);
void writeData (FILE *file, double t, double *y, int num_eq);
void plotSolution (int id);
