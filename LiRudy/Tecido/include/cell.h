#ifndef CELL_H_
#define CELL_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

// Condicoes iniciais

// Constantes
static constexpr double T_STIM = 0.5; 
static constexpr double PACING = 500;

// GEOMETRY
static constexpr double pi = 3.14;			
static constexpr double radius = 0.00175;	
static constexpr double length = 0.0164;	
static constexpr double rcg = 1.54;

class Cell
{
    
public:
    // State variables
    double v, m, h, j, d, f, f2, fca, fca2, xs1, xs2, xr, a, i, i2;
    double ml, ml3, hl, hl3, jl, jl3; 
    double casss, cajsr, cacsr, cansr, cassl;
    double nai, nassl, nasss, ki, cai, b, g, u, y, camktrap;

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

public:
    Cell ();
    void setCell ();
    void solve (int id, double t, double dt);
    void swap ();
    void write (double t, FILE *out);

    // TO DO
    friend ostream& operator<< (ostream &ost, const Cell &cell)
    {
        return ost;
    }

};

void computeGeometrics ();

void setBeats (int nbeats);

#endif