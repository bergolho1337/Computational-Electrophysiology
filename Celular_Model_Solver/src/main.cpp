#include <cstdlib>
#include <cstdio>
#include "../include/celularModel.h"

//define PLOT

using namespace std;

int main (int argc, char *argv[])
{
  printf("============================== CELULAR MODEL SOLVER =================================================\n");
  printf("[!] Solving the system of equations related to ionic currents (I_ion) using Explicit Euler\n");
  printf("-----------------------------------------------------------------------------------------------------\n");
  if (argc-1 < 3)
  {
    printf("Usage:> %s <dt> <t_max> <id_model>\n",argv[0]);
    printf("-----------------------------------------------------------------------------------------------------\n");
    printf("<dt> = Size of the discretization in time\n");
    printf("<t_max> = Time of the simulation\n");
    printf("<id_model> = Identification of the celular model\n");
    printf("\t1 - Mitchell & Shaeffer\n");
    printf("\t2 - Noble\n");
    printf("\t3 - Hodkin-Huxley\n");
    printf("\t4 - FitzHugh-Nagumo\n");
    printf("\t5 - Li & Rudy\n");
    printf("\t6 - Haq\n");
    printf("====================================================================================================\n");
    exit(-1);
  }
  else
  {
    CelularModel *cellModel = initModel(argc,argv);
    solveModel(cellModel);
    freeModel(cellModel);
    
    #ifdef PLOT
    plotSolution(cellModel->id);
    #endif
  }
}
