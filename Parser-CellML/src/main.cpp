#include <iostream>
#include "../include/parser.h"
#include "../include/timer.h"

int main (int argc, char *argv[])
{
  if (argc-1 < 1)
  {
    printf("=============== PARSER ========================\n");
    printf("Usage:> %s <input_file>\n",argv[0]);
    printf("<input_file> = Arquivo de entrada\n");
    printf("=============== PARSER ========================\n");
    exit(EXIT_FAILURE);
  }
  
  Parser *parser = new Parser(argc,argv);

  parser->convert();

  delete parser;
  
  return 0;
}
