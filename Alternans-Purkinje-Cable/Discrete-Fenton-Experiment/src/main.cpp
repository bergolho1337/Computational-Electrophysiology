#include <cstdio>
#include "../include/timer.h"
#include "../include/monodomain.h"

using namespace std;

int main (int argc, char *argv[])
{
  
  if (argc-1 != 1)
  {
    Usage(argv[0]);
    return 1;
  }
  else
  {
    double start, finish, elapsed;

    GET_TIME(start);
    solve_monodomain(argc,argv);
    GET_TIME(finish);
    elapsed = finish - start;
  
    printf("==========================================================\n");
    printf("[!] Time elapsed = %.10lf seconds\n",elapsed);
    printf("==========================================================\n");

    return 0;
  }
}