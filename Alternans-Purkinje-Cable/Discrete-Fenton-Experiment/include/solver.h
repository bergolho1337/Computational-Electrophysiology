#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
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

struct CVolume
{
  int type;                                   // 0 = Purkinje cell || 1 = PMJ
  double *yOld;                               // Solution at timestep n     {v^n, gate^n}
  double *yNew;                               // Solution at timestep n+1   {v^n+1, gate^n+1}
  double *yStar;                              // Solution at timestep n+1/2 {v^n+1/2, gate^n}
}typedef CVolume;

struct Derivative
{
  double t;                                   // Time of the maximum derivative
  double value;                               // Derivative value
}typedef Derivative;

struct Velocity
{
  FILE *velocityFile;                         // Reference to the file where the velocity will be stored
  int np;                                     // Number of volumes that the velocity will be calculated
  int id_source;                              // Identifier of the source volume
  int *ids;                                   // Identifier of the sink points (ids[0] -> source) 
  double t1;                                  // Time when the despolarization occur on the source volume 
  double *t2;                                 // Time when the despolarization occur on the sink volumes
}typedef Velocity;

// Structure for the plot volumes
struct Plot
{
  FILE **plotFile;                            // Reference to the file of the plot volume
  int np;                                     // Number of plot volumes
  int *ids;                                   // Identifier of the plot volume
}typedef Plot;

class Solver
{
    // Constants of the monodomain equation
    //static constexpr double BETA = 0.14;
    static constexpr double Cm = 1.2;
    //static constexpr double SIGMA = 0.004;
    static constexpr double h2 = 0.25;
    static constexpr double d2 = 0.002;
    static constexpr double RPMJ = 11000.0;
    static constexpr int OFFSET = 50;
public:
    Solver (User_Options *options);
    void solve ();
    void print ();
    void error (const char msg[]);
private:
    Graph *g;                           // Graph representing the Purkinje network
    CVolume *vol;                       // Vector of control volumes
    Derivative *dvdt;                   // Vector with the maximum derivative for each volume
    Plot *plot;                         // Vector with the plot ids
    Velocity *vel;                      // Vector with the propagation velocity of the plot ids
    Func *func;                         // Vector of function of the celullar model
    int M;                              // Number of timesteps in time
    double dt;                          // Size timestep in time
    double tmax;                        // Maximum time of the simulation
    double dx;                          // Size timestep in space
    string mesh_filename;               // Mesh filename
    string steady_filename;             // Input Steady-State filename
    string plot_filename;               // Plot id filename

    double alfa;                        // Parameter: R_pmj * Vol_pmj
    double d1;                          // Parameter: d1
    double BETA;                        // Surface / Volume ratio
    double SIGMA;                       // Conductivity citoplasm
    double GGAP;                        // Conductance of the gap junction

    void setSensibilityParam (User_Options *options);
    void setTypeCell ();  
    void setControlVolumes ();
    void setFunctions ();
    void setInitCondFromFile ();
    void setVelocityPoints ();
    void setDerivative ();
    void setPlot ();
    void setTerm ();
    void setMatrix (SpMat &a);
    void setMatrix2 (SpMat &a);
    void assembleLoadVector (VectorXd &b);
    void moveVstar (const VectorXd vm);
    bool isConnToPMJ (Edge *e);
    void solveODE (double t);
    void nextTimestep ();
    void calcMaxDerivative (double t);
    void calcVelocity ();
    void writeVTKFile (int iter);
    void writePlotData (double t);
    void writeStateVectorToFile (int count, int max_count, double t);
    void writeSteadyStateFile (FILE *sstFile);
};



void swap (double **a, double **b);
void printProgress (int iter, int max_iter);
int get_num_digits (int num);

#endif