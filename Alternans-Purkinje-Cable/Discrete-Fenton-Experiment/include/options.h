#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>

using namespace std;

class Options
{
public:
    Options (int argc, char *argv[]);
    void print_user_options ();
    
public:
    // Main section
    int num_threads;
    double dt;
    double tmax;
    bool steady_state;
    int print_rate;
    int sst_rate;
    string network_filename;
    string sst_filename;
    string plot_filename;

    // Cell section
    double start_h;
    int num_div_cell;
    double start_diameter;
    double sigma_c;
    double G_gap;

    // Stimulus section
    double stim_current;
    double stim_start;
    double stim_duration;
    double start_period;
    double end_period;
    double period_step;
    int n_cycles;
    int id_limit;
};

#endif