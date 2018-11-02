#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>

using namespace std;

class User_Options
{
public:
    User_Options (int argc, char *argv[]);
    void print_user_options ();
    
public:
    bool steady_state;
    double dt;
    double tmax;
    double start_h;
    int num_div_cell;
    string mesh_filename;
    string sst_filename;
    string plot_filename;
    double alfa;
    double diameter;
    double sigma_c;
    double G_gap;
};

#endif