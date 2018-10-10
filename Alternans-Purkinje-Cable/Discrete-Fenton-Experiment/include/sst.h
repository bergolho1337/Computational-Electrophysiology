#ifndef SST_H_
#define SST_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <omp.h>
#include <Eigen/Sparse>
#include "../include/purkinje.h"
#include "../include/noble.h"
#include "../include/options.h"

using namespace std;
using namespace Eigen;

// Macros for the Eigen library structures
typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;

// Function pointer
typedef double (*Func) (int type, int point, double t, double y[]);

struct Volume
{
  int type;                                   // 0 = Purkinje cell || 1 = PMJ
  double *yOld;                               // Solution at timestep n     {v^n, gate^n}
  double *yNew;                               // Solution at timestep n+1   {v^n+1, gate^n+1}
  double *yStar;                              // Solution at timestep n+1/2 {v^n+1/2, gate^n}
}typedef Volume;

class SteadyState
{
    // Constants of the monodomain equation
    //static constexpr double BETA = 0.14;
    static constexpr double Cm = 1.2;
    //static constexpr double SIGMA = 0.004;
    static constexpr double h2 = 0.25;
    static constexpr double d2 = 0.002;
    static constexpr double RPMJ = 11000.0;
public:
    SteadyState (User_Options *options);
    void solve ();
    void print ();
private:
    Graph *g;                           // Graph representing the Purkinje network
    Volume *vol;                        // Vector of control volumes
    Func *func;                         // Vector of function of the celullar model
    int M;                              // Number of timesteps in time
    double dt;                          // Size timestep in time
    double tmax;                        // Maximum time of the simulation
    double start_h;                     // Size of a cell
    double dx;                          // Discretization size
    int num_div_cell;                   // Number of divisions of a cell
    string mesh_filename;               // Mesh filename
    string steady_filename;             // Output Steady-State filename

    double alfa;                        // Parameter: R_pmj * Vol_pmj
    double d1;                          // Parameter: d1
    double BETA;                        // Surface / Volume ratio
    double SIGMA;                       // Conductivity of the citoplasm
    double GGAP;                        // Conductance of the gap junction

    void setSensibilityParam (User_Options *options);
    void setControlVolumes ();
    void setFunctions ();
    void setInitCondModel ();
    void setMatrix (SpMat &a);
    void setMatrix2 (SpMat &a);
    void assembleLoadVector (VectorXd &b);
    void moveVstar (const VectorXd vm);
    bool isConnToPMJ (Edge *e);
    void solveODE (double t);
    void nextTimestep ();
    void writeVTKFile (int iter);
    void writeSteadyStateFile (FILE *sstFile);

};

void swap (double **a, double **b);
void printProgress (int iter, int max_iter);

#endif