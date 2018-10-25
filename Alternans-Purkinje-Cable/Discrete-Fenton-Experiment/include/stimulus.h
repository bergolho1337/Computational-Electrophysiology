#ifndef STIMULUS_H_
#define STIMULUS_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "options.h"

using namespace std;

class Stimulus
{
public:
    Stimulus (Options *user_options);
    void print ();
public:
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