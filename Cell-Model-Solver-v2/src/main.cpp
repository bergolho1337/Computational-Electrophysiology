#include <cstdlib>
#include <cstdio>
#include "../include/celularModel.h"

using namespace std;

int main (int argc, char *argv[])
{
  printf("============================== CELULAR MODEL SOLVER =================================================\n");
  printf("[!] Solving the system of equations related to ionic currents (I_ion)\n");
  printf("-----------------------------------------------------------------------------------------------------\n");
  if (argc-1 != 4)
  {
    Usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  else
  {
    CelullarModel *cm = new CelullarModel(argc,argv);
    cm->Solve();
    delete cm;
  }
  return 0;
}
