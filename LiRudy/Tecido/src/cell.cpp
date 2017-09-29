#include "../include/cell.h"

// GEOMETRY
double ageo,acap,vcell,vmyo,vnsr,vjsr,vsss,vcsr,vssl,vmito;

Cell::Cell ()
{
    
}

// Compute the geometry related variables 
void computeGeometrics ()
{
    // CELL GEOMETRY
	vcell	= 1000*pi*radius*radius*length;
	ageo	= 2*pi*radius*radius + 2*pi*radius*length;
	acap	= rcg*ageo;
	vmyo	= vcell * 0.60;
	vnsr	= vcell * 0.04;
	vmito   = vcell * 0.18;
	vjsr	= vcell * 0.002;
	vcsr	= vcell * 0.008;
	vsss	= vcell * 0.02;
	vssl    = vcell * 0.15;
}

void Cell::setCell ()
{
    v		= -84.058830;
    m		= 0.000821;
    h		= 0.995741;
    j		= 0.999872;
    d		= 0.000016;
    f		= 0.999193;
    f2		= 0.988692;
    fca		= 0.965405;
    fca2	= 0.739378;
    xs1		= 0.001114;
    xs2		= 0.042234;
    xr		= 0.069808;
    a		= 0.000119;
    i		= 0.992541;
    i2		= 0.745628;
    ml		= 0.000329;
    ml3		= 0.046538;
    hl		= 0.984170;
    hl3		= 0.853893;
    jl		= 0.912569;
    jl3		= 0.827885;
    casss	= 0.000135;
    cajsr	= 1.510741;
    cacsr	= 1.537577;
    cansr	= 1.538668;
    cassl	= 0.000130; 
    nai		= 11.501546;
    nassl	= 11.501230;
    nasss	= 11.501240;
    ki		= 136.422946;
    cai		= 0.000053;
    b	    = 0.000437;
    g	    = 0.990384;
    u       = 0.535627;
    y       = 0.182859;
    camktrap= 0.010600;
}