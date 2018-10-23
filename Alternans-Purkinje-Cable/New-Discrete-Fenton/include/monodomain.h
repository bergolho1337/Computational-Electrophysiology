#ifndef MODEL_H_
#define MODEL_H_

#include <cstdio>
#include <cstdlib>
#include <string>

#include "options.h"
#include "solver.h"

class Monodomain
{
public:
    Monodomain (int argc, char *argv[]);
    void solve ();
    void error (const char msg[]);
private:
    Options *user_options;
    Solver *solver;
};

void Usage (const char p_name[]);
void solve_monodomain (int argc, char *argv[]);

#endif