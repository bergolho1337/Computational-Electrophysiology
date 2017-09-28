#include <iostream>
#include "../include/noble.h"
#include "../include/timer.h"

int main (int argc, char *argv[])
{
  if (argc-1 < 3)
  {
    printf("====== NOBLE 1962 ========\n");
    printf("Usage:> %s <t_max> <dt> <num_threads>\n",argv[0]);
    printf("<t_max> = Tempo maximo de simulacao (s)\n");
    printf("<dt> = Tamanho da discretizacao no tempo (s)\n");
    printf("<num_threads> = Numero de threads\n");
    exit(EXIT_FAILURE);
  }
  
  Noble *noble = new Noble(argc,argv);
  //cout << *noble << endl;
  double start, finish, elapsed;
  
  GET_TIME(start);
  noble->solve();
  GET_TIME(finish);
  elapsed = finish - start;
  cout << "Elapsed time = " << elapsed << " seconds" << endl;
  
  return 0;
}
