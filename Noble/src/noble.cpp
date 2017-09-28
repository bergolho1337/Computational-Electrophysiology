#include "../include/noble.h"

Noble::Noble (int argc, char *argv[])
{
  tmax = atof(argv[1]);
  dt = atof(argv[2]);
  nthreads = atoi(argv[3]);
  nbeats = tmax / PACING + 1;
  n = tmax / dt;

  allocMem();
  setInitCond();
  setBeats(nbeats);

}

void Noble::setInitCond ()
{
  for (int i = 0; i < (int)cells.size(); i++)
    cells[i].setCell();
}

void Noble::allocMem ()
{
  cells.assign(NCELLS,Cell(NEQ));
}

// Resolve o problema utilizando OpenMP
void Noble::solve ()
{
  double t;
  FILE *out = fopen("solution.dat","w+");

  for (int k = 0; k < n; k++)
  {
    t = k*dt;
    // Output
    cells[OUT_ID].write(t,out);
      
    // Cada thread fica responsavel por um conjunto de celulas
    #pragma omp parallel for num_threads(nthreads)  
    for (int i = 0; i < (int)cells.size(); i++)
      cells[i].solve(i,t,dt);

    for (int i = 0; i < (int)cells.size(); i++)
      cells[i].swap();
      
  }
  fclose(out);
}


/*
void readInput (int argc, char *argv[])
{
  t_max = atof(argv[1]);
  dt = atof(argv[2]);
}

void setInitialConditions (double *yOld)
{
  yOld[0] = v0;
  yOld[1] = m0;
  yOld[2] = h0;
  yOld[3] = n0;
  subindo = 0;
  descendo = 0;
}

void writeIteration (double t, double *y)
{
  FILE *file = fopen("solution.dat","a");
  int i;
  fprintf(file,"%e",t);
  for (i = 0; i < num_eq; i++)
    fprintf(file," %e",y[i]);
  fprintf(file,"\n");
  fclose(file);
}

void solveEDO ()
{
  int i, j, n;
  double t, f, *yOld, *yNew;

  // Numero de subintervalos
  n = nearbyint(t_max / dt);

  // Alocar vetor das solucoes
  yOld = (double*)malloc(sizeof(double)*num_eq);
  yNew = (double*)malloc(sizeof(double)*num_eq);

  // Setar as condicoes iniciais
  setInitialConditions(yOld);

  // Resolver a EDO
  for (i = 0; i < n; i++)
  {
    t = i*dt;
    // Escreve o valor da solucao no arquivo .dat
    writeIteration(t,yOld);
    for (j = 0; j < num_eq; j++)
    {
      // Calcular a funcao f
      switch(j)
      {
        case 0: f = dvdt(i,t,yOld);
                break;
        case 1: f = dmdt(t,yOld);
                break;
        case 2: f = dhdt(t,yOld);
                break;
        case 3: f = dndt(t,yOld);
                break;
      }
      // Metodo de Euler Explicito
      yNew[j] = yOld[j] + f*dt;
    }
    isEquilibrium(t,yOld,yNew);
    memcpy(yOld,yNew,sizeof(double)*num_eq);
  }

  // Libera memoria
  free(yOld);
  free(yNew);
}

void isEquilibrium (double t, double *yOld, double *yNew)
{
  if (yNew[0] > yOld[0])
  {
    if (descendo == 1)
      printf("[+] t = %e ... Despolarizacao de um PA :> v = %e || m = %e || h = %e || n = %e\n",t,yOld[0],yOld[1],yOld[2],yOld[3]);
    descendo = 0;
    subindo = 1;
  }
  else if (yNew[0] < yOld[0])
  {
    if (subindo == 1)
    {
      descendo = 1;
      subindo = 0;
    }
  }
}
*/