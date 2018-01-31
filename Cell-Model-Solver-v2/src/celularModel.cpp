#include "../include/celularModel.h"

CelullarModel::CelullarModel (int argc, char *argv[])
{
  DT = atof(argv[1]);
  TMAX = atof(argv[2]);
  ID = atoi(argv[3]);

  switch (ID)
  {
    case 1: cout << "[!] Building Mitchell and Shaeffer model ..." << endl;
            buildMitchell();
            break;
    case 2: cout << "[!] Building Noble model ..." << endl;
            buildNoble();
            break;
    case 3: cout << "[!] Building DiFrancesco model ..." << endl;
            buildDiFrancesco();
            break;
    
    default:  cout << "[-] Invalid ID code!" << endl;
              break;
  }
}

CelullarModel::~CelullarModel ()
{
  delete [] ALGEBRAIC;
  delete [] RATES;
  delete [] STATES;
  delete [] CONSTANTS;
}

void CelullarModel::buildMitchell ()
{
  NRATES = 2;
  NCONSTANTS = 10;
  NALGEBRAICS = 3;

  MODEL = new Mitchell("Mitchell");
  ALGEBRAIC = new double[NALGEBRAICS];
  RATES = new double[NRATES];
  STATES = new double[NRATES];
  CONSTANTS = new double[NCONSTANTS];

  MODEL->initConst(CONSTANTS,RATES,STATES);
}

void CelullarModel::buildNoble ()
{
  NRATES = 4;
  NCONSTANTS = 5;
  NALGEBRAICS = 12;

  MODEL = new Noble("Noble");
  ALGEBRAIC = new double[NALGEBRAICS];
  RATES = new double[NRATES];
  STATES = new double[NRATES];
  CONSTANTS = new double[NCONSTANTS];

  MODEL->initConst(CONSTANTS,RATES,STATES);
}

void CelullarModel::buildDiFrancesco ()
{
  NRATES = 16;
  NCONSTANTS = 50;
  NALGEBRAICS = 45;

  MODEL = new DiFrancesco("DiFrancesco");
  ALGEBRAIC = new double[NALGEBRAICS];
  RATES = new double[NRATES];
  STATES = new double[NRATES];
  CONSTANTS = new double[NCONSTANTS];

  MODEL->initConst(CONSTANTS,RATES,STATES);
}

void CelullarModel::Solve ()
{
  double start, end, elapsed;
  #ifdef OUTPUT
  ofstream out("Output/solution.dat");
  #endif

  int N = nearbyint(TMAX / DT);
  GET_TIME(start);
  for (int i = 0; i < N; i++)
  {
      double t = i*DT;
      MODEL->compVariables(t,CONSTANTS,RATES,STATES,ALGEBRAIC);
      MODEL->compRates(t,CONSTANTS,RATES,STATES,ALGEBRAIC);
      
      #ifdef OUTPUT
      PrintSolution(out,t,STATES);
      #endif

      for (int j = 0; j < NRATES; j++)
          STATES[j] = STATES[j] + RATES[j]*DT;
  }
  GET_TIME(end);
  elapsed = end - start;

  WriteLogFile(elapsed);
}

void CelullarModel::PrintSolution (ofstream &out, const double t, const double STATES[])
{
  out << fixed << setprecision(8) << t << " ";
  for (int i = 0; i < NRATES; i++)
      out << STATES[i] << " ";
  out << endl;
}

void Usage (const char pName[])
{
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "Usage:> " << pName << " <dt> <tmax> <id>" << endl;
  cout << "<dt> = Size of the discretization in time" << endl;
  cout << "<tmax> = Maximum simulation time" << endl;
  cout << "<id> = Identification of the celullar model" << endl;
  cout << "\t1 - Mitchell and Schaeffer" << endl;
  cout << "\t2 - Noble" << endl;
  cout << "\t3 - DiFrancesco" << endl;
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "Examples:> " << pName << " 0.1 500.0 1 (1 pulse Mitchell)" << endl;
  cout << "           " << pName << " 0.1 600 2 (1 pulse Noble)" << endl;
  cout << "           " << pName << " 1.0e-05 2.0 3 (1 pulse DiFrancesco)" << endl;
}

void CelullarModel::WriteLogFile (const double elapsedTime)
{
  ofstream logFile("Output/logfile.log");
  logFile << fixed << setprecision(6);
  switch (ID)
  {
    case 1:   logFile << "Mitchell" << endl;
              break;
    case 2:   logFile << "Noble" << endl;
              break;
    case 3:   logFile << "DiFrancesco" << endl;
              break;
    default:  cout << "[-] Invalid ID !" << endl;
              exit(EXIT_FAILURE);
  }
  logFile << elapsedTime << endl;
}

