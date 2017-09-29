#ifndef PRD_H
#define PRD_H

#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include "../include/cell.h"

using namespace std;

class LiRudy
{
    static constexpr int NCELLS = 10;               // Number of cells
    static constexpr int NSC = 2;                   // Number of stimulus cells 
    static constexpr int OUT_ID = 0;                // Plot cell identifier

private:
    int nthreads;
    int nbeats;
    int n;
    double dt;
    double tmax;
    vector<Cell> cells;

public:
    LiRudy (int argc, char *argv[]);
    void allocMem ();
    void setInitCond ();
    void setStimCells ();
    void solve ();

    // DEBUG
    friend ostream& operator<< (ostream &ost, const LiRudy &lirudy)
    {
        ost << "-------- PARAMETERS ------------" << endl;
        ost << "Dt = " << lirudy.dt << endl;
        ost << "tmax = " << lirudy.tmax << endl;
        ost << "Number of timesteps = " << lirudy.n << endl;
        ost << "--------------------------------" << endl;

        return ost;
    }

};

/* CONSTANTS */
static constexpr int MAX_FILENAME = 200;

// PHYSICAL CONSTANTS
static constexpr double frdy = 96485;		      
static constexpr double R = 8314;			  
static constexpr double temp = 310;		    

static constexpr double nao = 140;			
static constexpr double cao = 1.8;			
static constexpr double ko  = 5.4;			
static constexpr double clo = 100;		

static constexpr double zna = 1;			
static constexpr double zk  = 1;			
static constexpr double zcl = -1;		
static constexpr double zca = 2;			
static constexpr double ganai = 0.75;		
static constexpr double ganao = 0.75;		
static constexpr double gaki  = 0.75;	
static constexpr double gako  = 0.75;	
static constexpr double gacai = 1.0;		
static constexpr double gacao = 0.341;

// TIME
static constexpr double dt = 0.1;	
static constexpr double tmax = 5000;

// VOLTAGE			
static constexpr double dvdtthresh = 1;		

// STIMULUS CURRENT
static constexpr double stimdur = 0.5;

// MEMBRANE IONIC CURRENTS
static constexpr double gna = 18;	  
static constexpr double gnal2 = 0.052;	
static constexpr double gnal3 = 0.018;
static constexpr double pca = 1.9926e-4;	
static constexpr double powtau = 10;	
static constexpr double gcat = 0.07875;    	
static constexpr double gtos = 0.1414;	   
static constexpr double gtof = 0.042;  
static constexpr double prnak = 0.014;
static constexpr double gnab = 0.0025;     

static constexpr double pcab = 3.99e-8;	  
static constexpr double pnab = 0.64e-8;

static constexpr double inacamax = 2.52;
static constexpr double kmcaact = 0.000125;
static constexpr double kmnai1 = 12.3;		
static constexpr double kmnao = 87.5;		
static constexpr double kmcai = 0.0036;	
static constexpr double kmcao = 1.3;		
static constexpr double nu = 0.35;			
static constexpr double ksat = 0.27;	
static constexpr double ibarnak = 1.1004;
static constexpr double ipcabar = 0.0115;		   
static constexpr double kmpca = 0.0005;

// CALCIUM FLUXES RATE CONSTANTS
static constexpr double tautr1 = 120;
static constexpr double tautr2 = 120;	
static constexpr double gaptau = 12;
static constexpr double sstau = 0.2;

static constexpr double k1 = 150000;
static constexpr double k1a = 16.5;
static constexpr double k0 = 96000;
static constexpr double k0a = 9.6;
static constexpr double k2 = 1800;
static constexpr double k2a = 0.21;
static constexpr double tauip3r = 3.7;

static constexpr double dqupcamkbar = 0.75;
static constexpr double dkmplbbar = 0.00017;
static constexpr double nsrbar = 15.0;
static constexpr double bsrbar = 0.019975;	
static constexpr double kmbsr = 0.00087;		
static constexpr double bslbar = 0.4777;	
static constexpr double kmbsl = 0.0087;
static constexpr double csqnbar = 2.88;		    
static constexpr double kmcsqn = 0.8;

static constexpr double cmdnbar = 0.1125;	
static constexpr double kmcmdn = 2.38e-3;
static constexpr double trpnbar = 3.15e-2;
static constexpr double kmtrpn = 0.5e-3;

static constexpr double camk0 = 0.05;		
static constexpr double alphacamk = 0.05;		
static constexpr double betacamk = 0.00068;	
static constexpr double kmcam = 0.0015;		
static constexpr double kmcamk = 0.15;	
static constexpr double fca_dtaucamkbar = 10.0;

/*
struct Cell
{
    // State variables
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

    // Stimulus variables 
    int beats;
    int stimcount;
    double BCL;
    double S2;
    double tstim, stimtime, dvdtclock;

    // Currents
    double ina, inal, inal2, inal3, inab, inacass, inaca, inak, inatot;
    double ibarca;
    double icat, ical, icab, icatot;
    double itos, itof, ito1;
    double ikr, iks, ik1, iktot;
    double ipca;
    double ifna, ifk, iftotal;
    double itot;

    // Calcium dinamycs 
    double camkactive;
    double qip3, qrel1, qrel2, qup1, qup2, qtr1, qtr2;
    double caavg;

    // CAMKII DYNAMICS
    double camkbound;
    double fca_dtaucamk;

    // Calcium fluxes and concentrations
    double qdiff;
    double du,POip3;
    double irelss;
    double  dqupcamk;
    double  dkmplb;
    double bsss,csqn1,bjsr,cjsr,csqn,bcsr,ccsr;
    double dcasss,cassstot,bsr,bsl,b1,c1,d1;
    double dcassl;	
                    
    double dcajsr,cajsrtot;
    double dcacsr,cacsrtot;			         
    double dcansr;

    double dcai,catotal,cmdn;
    double trpn;
    double bmyo,cmyo,dmyo;				   


    // SODIUM/POTASSIUM FLUXES AND CONCENTRATIONS
    double dnai,dnasss,dki,ksss,dksss;	
    double qgap,qdiffna,qgapna,dnassl;

    // Reverse potential
    double ena, ek, eca;

    // Membrane ionic currents
    double ma,mb,mtau,mss,ha,hb,htau,hss,ja,jb,jtau,jss;
    double alphaml,betaml,mltau,mlss,hltau,hlss;
    double i3tau,i3ss,i3,Rtau,Rss,Ri,ml3tau,ml3ss,hl3tau,hl3ss,ireltau,REL;
    double jltau,jlss,jl3tau,jl3ss;
    double dss,dtau,dpower,powss;
    double fss,ftau,f2ss,f2tau,fcass,fcatau,fca2tau,fca2ss;
    double taub,taug,bss,gss;
    double rto1,alphaa,betaa,atau,ass,alphai,betai,itau,iss,alphai2,betai2,i2tau,i2ss;
    double gkr,xrss,xrtau,rkr;
    double gks,eks;
    double xs1tau,xs2tau,xsss;
    double gk1,k1ss;
    double yss, ytau;
    double allo,num,denommult,denomterm1,denomterm2,deltaE;

    double dvdt;


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
*/

#endif