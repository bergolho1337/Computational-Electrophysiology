#include "../include/celularModel.h"

CelullarModel::CelullarModel (int argc, char *argv[])
{
  DT = atof(argv[1]);
  TMAX = atof(argv[2]);
  ID = atoi(argv[3]);
  SOLVER = atoi(argv[4]);

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
    case 4: cout << "[!] Building LuoRudy model ..." << endl;
            buildLuoRudy();
            break;
    case 5: cout << "[!] Building Fitz-Hugh Nagumo model ..." << endl;
            buildFitzHugh();
            break;
    case 6: cout << "[!] Building TenTusscher Miocardium Cell model ..." << endl;
            buildTenTusscherMC();
            break;
    case 7: cout << "[!] Building TenTusscher Epicardium Cell model ..." << endl;
            buildTenTusscherEPI();
            break;
    case 8: cout << "[!] Building TenTusscher Endocardium Cell model ..." << endl;
            buildTenTusscherENDO();
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

void CelullarModel::buildLuoRudy ()
{
  NRATES = 8;
  NCONSTANTS = 24;
  NALGEBRAICS = 25;

  MODEL = new LuoRudy("LuoRudy");
  ALGEBRAIC = new double[NALGEBRAICS];
  RATES = new double[NRATES];
  STATES = new double[NRATES];
  CONSTANTS = new double[NCONSTANTS];

  MODEL->initConst(CONSTANTS,RATES,STATES);
}

void CelullarModel::buildFitzHugh ()
{
  NRATES = 2;
  NCONSTANTS = 3;
  NALGEBRAICS = 1;

  MODEL = new FitzHugh("Fitz-Hugh");
  ALGEBRAIC = new double[NALGEBRAICS];
  RATES = new double[NRATES];
  STATES = new double[NRATES];
  CONSTANTS = new double[NCONSTANTS];

  MODEL->initConst(CONSTANTS,RATES,STATES);
}

void CelullarModel::buildTenTusscherMC ()
{
  NRATES = 19;
  NCONSTANTS = 53;
  NALGEBRAICS = 70;

  MODEL = new TenTusscherMC("TenTusscherMC");
  ALGEBRAIC = new double[NALGEBRAICS];
  RATES = new double[NRATES];
  STATES = new double[NRATES];
  CONSTANTS = new double[NCONSTANTS];

  MODEL->initConst(CONSTANTS,RATES,STATES);
}

void CelullarModel::buildTenTusscherEPI ()
{
  NRATES = 19;
  NCONSTANTS = 53;
  NALGEBRAICS = 70;

  MODEL = new TenTusscherEPI("TenTusscherEPI");
  ALGEBRAIC = new double[NALGEBRAICS];
  RATES = new double[NRATES];
  STATES = new double[NRATES];
  CONSTANTS = new double[NCONSTANTS];

  MODEL->initConst(CONSTANTS,RATES,STATES);
}

void CelullarModel::buildTenTusscherENDO ()
{
  NRATES = 19;
  NCONSTANTS = 53;
  NALGEBRAICS = 70;

  MODEL = new TenTusscherENDO("TenTusscherENDO");
  ALGEBRAIC = new double[NALGEBRAICS];
  RATES = new double[NRATES];
  STATES = new double[NRATES];
  CONSTANTS = new double[NCONSTANTS];

  MODEL->initConst(CONSTANTS,RATES,STATES);
}

void CelullarModel::RushLarsen ()
{
  // Mitchell-Schaeffer and Fitz-Hugh Nagumo do not have the gating variables
  if (ID == 1 || ID == 5)
  {
    cerr << "[-] ERROR ! The model cannot be solved with Rush-Larsen" << endl;
    cerr << "[*] Try solving the model using the Forward Euler integrator" << endl;
    exit(EXIT_FAILURE);
  }

  cout << "[!] Solving celular model using Rush-Larsen !" << endl;
  
  double start, end, elapsed;
  #ifdef OUTPUT
  ofstream out("Output/solution.dat");
  #endif

  int N = nearbyint(TMAX / DT);
  GET_TIME(start);
  for (int i = 0; i < N; i++)
  {
      double t = i*DT;
      MODEL->compRates_FE(t,CONSTANTS,RATES,STATES,ALGEBRAIC);
      MODEL->compVariables(t,CONSTANTS,RATES,STATES,ALGEBRAIC);
      
      #ifdef OUTPUT
      PrintSolution(out,i,STATES);
      #endif

      // Calls how to solve the Rush-Larsen for the celular model
      MODEL->RL(DT,CONSTANTS,RATES,STATES,ALGEBRAIC);
      
  }
  GET_TIME(end);
  elapsed = end - start;

  WriteLogFile(elapsed);
}

void CelullarModel::ForwardEuler ()
{
  
  cout << "[!] Solving celular model using Forward Euler !" << endl;

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
      PrintSolution(out,i,STATES);
      #endif

      for (int j = 0; j < NRATES; j++)
          STATES[j] = STATES[j] + RATES[j]*DT;
  }
  GET_TIME(end);
  elapsed = end - start;

  WriteLogFile(elapsed);
}

void CelullarModel::Solve ()
{
  switch (SOLVER)
  {
    case 1: ForwardEuler();
            break;
    case 2: RushLarsen();
            break;
  }
}

void CelullarModel::PrintSolution (ofstream &out, const int k, const double STATES[])
{
  double t = k*DT;
  if (k % PRINTRATE == 0)
  {
    out << fixed << setprecision(8) << t << " ";
    for (int i = 0; i < NRATES; i++)
      out << STATES[i] << " ";
    out << endl;
  }
}

void Usage (const char pName[])
{
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "Usage:> " << pName << " <dt> <tmax> <id_cell_model> <id_solver>" << endl;
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "<dt> = Size of the discretization in time" << endl;
  cout << "<tmax> = Maximum simulation time" << endl;
  cout << "<id_cell_model> = Identification of the celullar model" << endl;
  cout << "\t1 - Mitchell and Schaeffer" << endl;
  cout << "\t2 - Noble" << endl;
  cout << "\t3 - DiFrancesco" << endl;
  cout << "\t4 - Luo and Rudy 1" << endl;
  cout << "\t5 - Fitz-Hugh Nagumo" << endl;
  cout << "\t6 - TenTusscher Miocardium Cell" << endl;
  cout << "\t7 - TenTusscher Epicardium Cell" << endl;
  cout << "\t8 - TenTusscher Endocardium Cell" << endl;
  cout << "<id_solver> = Identification of the integrator" << endl;
  cout << "\t1 - Forward Euler" << endl;
  cout << "\t2 - Rush-Larsen" << endl;
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "Examples:> " << pName << " 0.1 500.0 1 1 (1 pulse Mitchell)" << endl;
  cout << "           " << pName << " 0.1 600 2 1 (1 pulse Noble)" << endl;
  cout << "           " << pName << " 0.01 2000 3 1 (1 pulse DiFrancesco)" << endl;
  cout << "           " << pName << " 0.01 2000 4 1 (1 pulse LuoRudy)" << endl;
  cout << "           " << pName << " 0.01 1000 5 1 (5 pulses Fitz-Hugh)" << endl;
  cout << "           " << pName << " 0.01 500 6 2 (1 pulse TenTusscher MCell)" << endl;
  cout << "           " << pName << " 0.01 500 7 2 (1 pulse TenTusscher EpiCell)" << endl;
  cout << "           " << pName << " 0.01 500 8 2 (1 pulse TenTusscher EndoCell)" << endl;
  
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
    case 4:   logFile << "LuoRudy" << endl;
              break;
    case 5:   logFile << "FitzHugh" << endl;
              break;
    case 6:   logFile << "TenTusscherMC" << endl;
              break;
    case 7:   logFile << "TenTusscherEPI" << endl;
              break;
    case 8:   logFile << "TenTusscherENDO" << endl;
              break;
    default:  cout << "[-] Invalid ID !" << endl;
              exit(EXIT_FAILURE);
  }
  logFile << elapsedTime << endl;
}

