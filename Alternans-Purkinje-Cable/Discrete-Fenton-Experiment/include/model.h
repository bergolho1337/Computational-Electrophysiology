#ifndef MODEL_H_
#define MODEL_H_

#include <cstdio>
#include <cstdlib>
#include <string>
#include "../include/sst.h"
#include "../include/solver.h"
#include "../include/options.h"
#include "../include/purkinje.h"

using namespace std;

class Model
{
public:
    Model (int argc, char *argv[]);
    void solve ();
    void error (const char msg[]);
private:
    SteadyState *sst;
    Solver *sol;
};

void Usage (const char pName[]);
void solveModel (int argc, char *argv[]);

#endif