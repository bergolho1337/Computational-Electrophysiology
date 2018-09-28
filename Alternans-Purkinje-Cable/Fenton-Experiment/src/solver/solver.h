#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <omp.h>
#include <getopt.h>
#include "../ini_parser/ini_file_sections.h"
#include "../ini_parser/ini.h"
#include "../stimuli/stimuli.h"
#include "../purkinje/purkinje.h"

//#include <Eigen/Sparse>
//#include "../graph/graph.h"
//#include "../noble/noble.h"

using namespace std;
//using namespace Eigen;

// Macros for the Eigen library structures
//typedef SparseMatrix<double> SpMat;
//typedef Triplet<double> T;

// Function pointer
//typedef double (*Func) (int type, int point, double t, double y[]);

/*
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
    static constexpr int OFFSET = 5;
public:
    Solver (int argc, char *argv[]);
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
    double SIGMA;                       // Conductivity Gap + Citoplasm

    void setSensibilityParam (int argc, char *argv[]);
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
};
*/

#define PRINT_LINE "-------------------------------------------------------------------------------------"
#define PRINT_LINE_2 "====================================================================================="

#define SIGMA 1400      // Dummy values
#define STIM_OPT 2000   // Dummy values
#define BETA 4000       // Dummy values
#define CM 5000         // Dummy values

struct user_options
{
  int num_threads;
  bool num_threads_was_set;
  double final_time;
  bool final_time_was_set;
  double dt;
  bool dt_was_set;
  int print_rate;
  bool print_rate_was_set;
  char *out_dir_name;
  bool out_dir_name_was_set;
  char *out_steady_state_dir;
  bool out_steady_state_dir_was_set;
  char *steady_state_filename;
  bool steady_state_filename_was_set;
  double sigma;
  bool sigma_was_set;
  double beta;
  bool beta_was_set;
  double cm;
  bool cm_was_set;

  char *config_file;

  struct stim_config *stim_configs;
  struct purkinje_config *purkinje_config;
};

struct user_options* new_user_options (int argc, char *argv[]);

void usage (const char pname[]);
void solve_model (int argc, char *argv[]);
void get_config_file (int argc, char *argv[], struct user_options *user_args);
int parse_config_file (void *user, const char *section, const char *name, const char *value);
void swap (double **a, double **b);
void printProgress (int iter, int max_iter);


#endif