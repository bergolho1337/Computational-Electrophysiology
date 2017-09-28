#ifndef CELL_H_
#define CELL_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

// Condicoes iniciais
//static constexpr double v0 = -8.0e+01;				 // mV
//static constexpr double m0 = 1.0e-02;
//static constexpr double h0 = 8.0e-01;
//static constexpr double n0 = 1.0e-02;

static constexpr double v0 = -8.290667e+01;
static constexpr double m0 = 4.029594e-02;
static constexpr double h0 = 8.783757e-01;
static constexpr double n0 = 5.487358e-01;

// Constantes
static constexpr double T_STIM = 2.0; 
static constexpr double V_STIM = 500;
static constexpr double V_GATE = -87;
static constexpr double PACING = 500;
static constexpr double G_NA_MAX = 4.0e+02;
static constexpr double E_NA = 4.0e+01;
static constexpr double G_L = 7.5e-02;
static constexpr double E_L = -6.0e+01;
static constexpr double CM = 1.2e+01;

class Cell
{
    
public:
    vector<double> yOld;
    vector<double> yNew;
public:
    Cell (int neq);
    void setCell ();
    void solve (int id, double t, double dt);
    void swap ();
    void write (double t, FILE *out);

    friend ostream& operator<< (ostream &ost, const Cell &cell)
    {
        for (int i = 0; i < (int)cell.yOld.size()-1; i++)
            ost << cell.yOld[i] << " ";
        ost << cell.yOld.back() << endl;
        return ost;
    }

    // dV/dt
    double I_Stim (int k, double t);
    double I_Leak (double t, double *y);
    double g_K2 (double t, double *y);
    double g_K1 (double t, double *y);
    double I_K (double t, double *y);
    double g_Na (double t, double *y);
    double I_Na (double t, double *y);
    double dvdt (int k, double t, double *y);
    
    // dm/dt
    double beta_m (double t, double *y);
    double alpha_m (double t, double *y);
    double dmdt (int k, double t, double *y);
    
    // dh/dt
    double beta_h (double t, double *y);
    double alpha_h (double t, double *y);
    double dhdt (int k, double t, double *y);
    
    // dn/dt
    double beta_n (double t, double *y);
    double alpha_n (double t, double *y);
    double dndt (int k, double t, double *y);
};

void setBeats (int nbeats);

#endif