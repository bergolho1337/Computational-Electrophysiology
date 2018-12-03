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

#include "constants.h"
#include "options.h"
#include "stimulus.h"
#include "purkinje.h"
#include "noble.h"

//#include "../include/purkinje.h"
//#include "../include/noble.h"
//#include "../include/options.h"

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
    //static constexpr double Cm = 1.2;
    //static constexpr double SIGMA = 0.004;
    //static constexpr double h2 = 0.25;
    //static constexpr double d2 = 0.002;
    //static constexpr double RPMJ = 11000.0;
    //static constexpr int OFFSET = 50;
public:
    Solver (Options *user_options);
    void solve ();
    void print ();
    //void error (const char msg[]);
private:
    int M;                              // Number of timesteps in time
    double dt;                          // Size timestep in time
    double tmax;                        // Maximum time of the simulation
    double dx;                          // Size timestep in space
    string network_filename;            // Network filename
    string steady_filename;             // Input Steady-State filename
    string plot_filename;               // Plot id filename
    int nthreads;                       // Number of threads
    int print_rate;                     // Rate which the VTK files will be saved
    int sst_rate;                       // Rate which the SST file will be saved

    double diameter;                    // Diameter of the Purkinje cell
    double beta;                        // Surface per Volume ratio
    double sigma_c;                     // Conductivity citoplasm
    double g_GAP;                       // Conductance of the gap junction
    
    Graph *the_purkinje_network;        // Graph representing the Purkinje network
    CVolume *vol;                       // Vector of control volumes
    Func *func;                         // Vector of function of the celular model
    Derivative *dvdt;                   // Vector with the maximum derivative for each volume
    Plot *plot;                         // Vector with the plot ids
    Velocity *vel;                      // Vector with the propagation velocity of the plot ids
    Stimulus *stim_config;              // Stimulus configuration

    void set_plot_cells ();
    void set_control_volumes ();
    void set_derivative ();
    void set_functions ();
    void set_velocity_points ();
    void set_plot_points ();
    void set_initial_conditions_from_file ();
    void set_initial_conditions_from_default ();

    void set_matrix (SpMat &A);
    void assemble_load_vector (VectorXd &b);
    void move_Vstar (const VectorXd vm);
    void solve_ODE (double t);
    void calc_max_derivative (double t, double current_period);
    void calc_velocity ();
    void next_timestep ();

    void write_plot_data (double t);
    void write_VTK_file (int iter);
    void write_steady_state_file (FILE *sstFile);

    /*
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
    */
};

void swap (double **a, double **b);
void print_progress (int iter, int max_iter);
int get_num_digits (int num);

#endif