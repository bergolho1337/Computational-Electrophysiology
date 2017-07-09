#ifndef PRD_H
#define PRD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* CONSTANTS */
const int MAX_FILENAME = 200;

// GEOMETRY
const double pi = 3.14;			
const double radius = 0.00175;	
const double length = 0.0164;	
const double rcg = 1.54;

// PHYSICAL CONSTANTS
const double frdy = 96485;		      
const double R = 8314;			  
const double temp = 310;		    

const double nao = 140;			
const double cao = 1.8;			
const double ko  = 5.4;			
const double clo = 100;		

const double zna = 1;			
const double zk  = 1;			
const double zcl = -1;		
const double zca = 2;			
const double ganai = 0.75;		
const double ganao = 0.75;		
const double gaki  = 0.75;	
const double gako  = 0.75;	
const double gacai = 1.0;		
const double gacao = 0.341;

// TIME
const double dt = 0.1;	
const double tmax = 5000;

// VOLTAGE			
const double dvdtthresh = 1;		

// STIMULUS CURRENT
const double stimdur = 0.5;

// MEMBRANE IONIC CURRENTS
const double gna = 18;	  
const double gnal2 = 0.052;	
const double gnal3 = 0.018;
const double pca = 1.9926e-4;	
const double powtau = 10;	
const double gcat = 0.07875;    	
const double gtos = 0.1414;	   
const double gtof = 0.042;  
const double prnak = 0.014;
const double gnab = 0.0025;     

const double pcab = 3.99e-8;	  
const double pnab = 0.64e-8;

const double inacamax = 2.52;
const double kmcaact = 0.000125;
const double kmnai1 = 12.3;		
const double kmnao = 87.5;		
const double kmcai = 0.0036;	
const double kmcao = 1.3;		
const double nu = 0.35;			
const double ksat = 0.27;	
const double ibarnak = 1.1004;
const double ipcabar = 0.0115;		   
const double kmpca = 0.0005;

// CALCIUM FLUXES RATE CONSTANTS
const double tautr1 = 120;
const double tautr2 = 120;	
const double gaptau = 12;
const double sstau = 0.2;

const double k1 = 150000;
const double k1a = 16.5;
const double k0 = 96000;
const double k0a = 9.6;
const double k2 = 1800;
const double k2a = 0.21;
const double tauip3r = 3.7;

const double  dqupcamkbar = 0.75;
const double  dkmplbbar = 0.00017;
const double nsrbar = 15.0;
const double bsrbar = 0.019975;	
const double kmbsr = 0.00087;		
const double bslbar = 0.4777;	
const double kmbsl = 0.0087;
const double csqnbar = 2.88;		    
const double kmcsqn = 0.8;

const double cmdnbar = 0.1125;	
const double kmcmdn = 2.38e-3;
const double trpnbar = 3.15e-2;
const double kmtrpn = 0.5e-3;

const double camk0 = 0.05;		
const double alphacamk = 0.05;		
const double betacamk = 0.00068;	
const double kmcam = 0.0015;		
const double kmcamk = 0.15;	
const double fca_dtaucamkbar = 10.0;

struct Cell
{
    /* State variables */
    double v;
    double m;
    double h;
    double j;
    double d;
    double f;
    double f2;
    double fca;
    double fca2;
    double xs1;
    double xs2;
    double xr;
    double a;
    double i;
    double i2;
    double ml;
    double ml3;
    double hl;
    double hl3;
    double jl;
    double jl3;
    double casss;
    double cajsr;
    double cacsr;
    double cansr;
    double cassl; 
    double nai;
    double nassl;
    double nasss;
    double ki;
    double cai;
    double b;
    double g;
    double u;
    double y;
    double camktrap;

    /* Stimulus variables */
    int beats;
    int stimcount;
    double BCL;
    double S2;
    double tstim, stimtime, dvdtclock;

    /* Currents */
    double ina, inal, inal2, inal3, inab, inacass, inaca, inak, inatot;
    double icat, ical, icab, icatot;
    double itos, itof, ito1;
    double ikr, iks, ik1, iktot;
    double ipca;
    double ifna, ifk, iftotal;
    double itot;

    /* Calcium dinamycs */
    double camkactive;
    double qip3, qrel1, qrel2, qup1, qup2, qtr1, qtr2;
    double caavg;


}typedef Cell;

FILE** createFiles (int n);
void setTimeSettings (Cell c[], int n);
void setInitialConditions (Cell c[], int n);
void setStimulusCells (Cell c[], int nsc, int n);
void compGeometrics ();

void solveModel (Cell c[], int n);

void timestep (Cell *c, int n);

void comp_revs (Cell *c);
void comp_ina (Cell *c);
void comp_inal (Cell *c);
void comp_inab (Cell *c);
void comp_ical (Cell *c);
void comp_icat (Cell *c);
void comp_icab (Cell *c);
void comp_ito1 (Cell *c);
void comp_ikr (Cell *c);
void comp_iks (Cell *c);
void comp_ik1 (Cell *c);
void comp_inaca (Cell *c);
void comp_inak (Cell *c);
void comp_ipca (Cell *c);
void comp_if (Cell *c);
void comp_istim (Cell *c);
void comp_itot (Cell *c);

void comp_ip3 (Cell *c);
void comp_qrel1 (Cell *c);
void comp_qrel2 (Cell *c);
void comp_qup1 (Cell *c);
void comp_qup2 (Cell *c);
void comp_qtr1 (Cell *c);
void comp_qtr2 (Cell *c);
void comp_conc (Cell *c);

void printtofile (FILE *ap, Cell c);

void printCells (Cell c[], int n);

#endif