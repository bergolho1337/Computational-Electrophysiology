#include "celularModel.h"

CelularModel* initModel (int argc, char *argv[])
{
  CelularModel *cellModel = (CelularModel*)malloc(sizeof(CelularModel));
  cellModel->dt = atof(argv[1]);
  cellModel->t_max = atof(argv[2]);
  cellModel->id = atoi(argv[3]);

  switch(cellModel->id)
  {
    case 1: printf("[!] Building Mitchell & Shaeffer model ... ");
            fflush(stdin);
            buildMitchell(cellModel);
            break;
    case 2: printf("[!] Building Noble model ... ");
            fflush(stdin);
            buildNoble(cellModel);
            break;
    case 3: printf("[!] Building Hodkin-Huxley model ... ");
            fflush(stdin);
            buildHodkin(cellModel);
            break;
    case 4: printf("[!] Building FitzHugh-Nagumo model ... ");
            fflush(stdin);
            buildFitzhugh(cellModel);
            break;
    case 5: printf("[!] Building Li & Rudy model ... ");
            fflush(stdin);
            buildLiRudy(cellModel);
            break;
    case 6: printf("[!] Building Haq model ... ");
            fflush(stdin);
            buildHaq(cellModel);
            break;
  }

  printf("ok\n");
  return cellModel;
}

void buildMitchell (CelularModel *cm)
{
  cm->num_eq = 2;
  cm->y0 = new double[cm->num_eq];
  cm->f = new Func[cm->num_eq];
  setInitialConditions__Mit(cm->y0,cm->num_eq);
  setFunctions__Mit(cm->f,cm->num_eq);
}

void buildNoble (CelularModel *cm)
{
  cm->num_eq = 4;
  cm->y0 = new double[cm->num_eq];
  cm->f = new Func[cm->num_eq];
  setInitialConditions__Nob(cm->y0,cm->num_eq);
  setFunctions__Nob(cm->f,cm->num_eq);
}

// Ta bugado !!!!
void buildHodkin (CelularModel *cm)
{
  cm->num_eq = 4;
  cm->y0 = new double[cm->num_eq];
  cm->f = new Func[cm->num_eq];
  setInitialConditions__Hod(cm->y0,cm->num_eq);
  setFunctions__Hod(cm->f,cm->num_eq);
}

void buildFitzhugh (CelularModel *cm)
{
  cm->num_eq = 2;
  cm->y0 = new double[cm->num_eq];
  cm->f = new Func[cm->num_eq];
  setInitialConditions__Fitz(cm->y0,cm->num_eq);
  setFunctions__Fitz(cm->f,cm->num_eq);
}

void buildLiRudy (CelularModel *cm)
{

}

void buildHaq (CelularModel *cm)
{

}

void plotSolution (int id)
{
  char cmd[50];
  int ret;
  sprintf(cmd,"python plot.py %d",id);
  ret = system(cmd);
  if (ret)
    printf("[-] ERROR! Plotting solution.\n");
}

void writeData (FILE *file, double t, double *y, int num_eq)
{
  fprintf(file,"%e",t);
  for (int i = 0; i < num_eq; i++)
    fprintf(file," %e",y[i]);
  fprintf(file,"\n");
}

void solveModel (CelularModel *cm)
{
  FILE *file;
  int n;
  double t, eval;
  double *yOld, *yNew;

  file = fopen("data.dat","w+");
  n = nearbyint(cm->t_max / cm->dt);
  yOld = new double[cm->num_eq];
  yNew = new double[cm->num_eq];
  memcpy(yOld,cm->y0,sizeof(double)*cm->num_eq);

  for (int k = 0; k < n; k++)
  {
    t = k*cm->dt;
    writeData(file,t,yOld,cm->num_eq);
    for (int i = 0; i < cm->num_eq; i++)
    {
      eval = cm->f[i](k,t,yOld);
      yNew[i] = yOld[i] + eval*cm->dt;
    }
    memcpy(yOld,yNew,sizeof(double)*cm->num_eq);
  }
  fclose(file);
  delete [] yOld;
  delete [] yNew;
}

void freeModel (CelularModel *cm)
{
  printf("[!] Desallocating celular model ... ");
  fflush(stdin);
  delete [] cm->y0;
  delete [] cm->f;
  delete cm;
  printf("ok\n");
}
