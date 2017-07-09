#include "prd.h"

// GEOMETRY
double ageo,acap,vcell,vmyo,vnsr,vjsr,vsss,vcsr,vssl,vmito;

// REVERSAL POTENTIALS
double ena,ek,eca;

// TIME
double t, T;

// VOLTAGE
double v,dvdt,dvdtclock;

// STIMULUS CURRENT
double istim = -80;	        
double tstim,stimtime;			   
int stimcount;

// TOTAL TRANSMEMBRANE CURRENTS
double icatot,iktot,inatot,icltot,itot;

// MEMBRANE IONIC CURRENTS
double ina;
double m,ma,mb,mtau,mss,h,ha,hb,htau,hss,j,ja,jb,jtau,jss;
double inal;
double ml,alphaml,betaml,mltau,mlss,hl,hltau,hlss;	
double jltau, jlss, jl, jl3tau, jl3ss, jl3;

double inal2, inal3;
double i3tau, i3ss, i3,Rtau, Rss, Ri,ml3tau, ml3ss, hl3tau, hl3, ml3, hl3ss,ireltau, REL;			

double ical,ibarca;
double d,dss,dtau,dpower,powss;
double f,fss,ftau,f2,f2ss,f2tau,fca,fcass,fcatau,fca2,fca2tau,fca2ss;

double icat;
double taub,taug,bss,gss,b,g;

double icab;

double ito1;			
double itos, itof;
double rto1,a,alphaa,betaa,atau,ass,i,alphai,betai,itau,iss,i2,alphai2,betai2,i2tau,i2ss;

double ikr,gkr,xr,xrss,xrtau,rkr;

double iks,gks,eks;
double xs1,xs1tau,xs2,xs2tau,xsss;

double ik1,gk1,k1ss;
double yss, ytau, y, ifna, ifk, iftotal;

double inab;

double inaca;	
double allo,num,denommult,denomterm1,denomterm2,deltaE,inacass;

double inak;

double ipca;

// CALCIUM FLUXES AND CONCENTRATIONS
double qrel1,qrel2;			
double irelss;	
double IP3 = 0.0001;
double du,u,POip3,qip3;

double  dqupcamk;
double  dkmplb,qup1,qup2;
double kmup   = 0.00028;

double  qtr1,qtr2;

double bsss,csqn1,bjsr,cjsr,csqn,bcsr,ccsr;
double casss,dcasss,cassstot,bsr,bsl,b1,c1,d1;
double cassl,dcassl;	
double qdiff;				               

double cajsr,dcajsr,cajsrtot;
double cacsr,dcacsr,cacsrtot;			         
double cansr,dcansr;

double cai,dcai,catotal,cmdn;
double trpn;
double bmyo,cmyo,dmyo;				   
double trpnbar1 = 3.5e-3;
double cmdnbar1 = 1.25e-2;
double csqnbar1 = 1.2;

double caavg;

// SODIUM/POTASSIUM FLUXES AND CONCENTRATIONS
double nai,dnai,nasss,dnasss,ki,dki,ksss,dksss;	
double qgap, qdiffna, qgapna, nassl, dnassl;

// CAMKII DYNAMICS
double camkbound,camkactive,camktrap;
double fca_dtaucamk;

// OUTPUT FILE
int count;

/* Compute the geometry related variables */
void compGeometrics ()
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

void setTimeSettings (Cell c[], int n)
{
    int i;
    t = 0;
	T = 0;
    for (i = 0; i < n; i++)
    {
        c[i].dvdtclock = 1000;
        c[i].stimtime = 1000;
        c[i].stimcount 	= -1;
    }
}

void setInitialConditions (Cell c[], int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        c[i].v		= -84.058830;
        c[i].m		= 0.000821;
        c[i].h		= 0.995741;
        c[i].j		= 0.999872;
        c[i].d		= 0.000016;
        c[i].f		= 0.999193;
        c[i].f2		= 0.988692;
        c[i].fca		= 0.965405;
        c[i].fca2	= 0.739378;
        c[i].xs1		= 0.001114;
        c[i].xs2		= 0.042234;
        c[i].xr		= 0.069808;
        c[i].a		= 0.000119;
        c[i].i		= 0.992541;
        c[i].i2		= 0.745628;
        c[i].ml		= 0.000329;
        c[i].ml3		= 0.046538;
        c[i].hl		= 0.984170;
        c[i].hl3		= 0.853893;
        c[i].jl		= 0.912569;
        c[i].jl3		= 0.827885;
        c[i].casss	= 0.000135;
        c[i].cajsr	= 1.510741;
        c[i].cacsr	= 1.537577;
        c[i].cansr	= 1.538668;
        c[i].cassl	= 0.000130; 
        c[i].nai		= 11.501546;
        c[i].nassl	= 11.501230;
        c[i].nasss	= 11.501240;
        c[i].ki		= 136.422946;
        c[i].cai		= 0.000053;
        c[i].b	    = 0.000437;
        c[i].g	    = 0.990384;
        c[i].u       = 0.535627;
        c[i].y       = 0.182859;
        c[i].camktrap= 0.010600;
    }
}

void setStimulusCells (Cell c[], int nsc, int n)
{
    int i;
    for (i = 0; i < nsc; i++)
    {
        c[i].beats = 10;
        c[i].BCL = 1000;
        c[i].S2 = 1000;
        c[i].tstim = 0;
    }
    for (i = nsc; i < n; i++)
    {
        c[i].beats = 10;
        c[i].BCL = tmax + 1000;
        c[i].S2 = 1000;
        c[i].tstim = c[i].BCL;
    }
}

void solveModel (Cell c[], int n)
{
    FILE **files = createFiles(n);
    int i;
    // Time loop
    while (t <= tmax)
    {
        timestep(c,n);
        // Space loop
        for (i = 0; i < n; i++)
        {
            // Compute each current of the model
            comp_revs(&c[i]);
            comp_ina (&c[i]);
            comp_inal (&c[i]);
            comp_inab (&c[i]);
            comp_ical (&c[i]);
            comp_icat (&c[i]);
            comp_icab (&c[i]);
            comp_ito1 (&c[i]);
            comp_ikr (&c[i]);
            comp_iks (&c[i]);
            comp_ik1 (&c[i]);
            comp_inaca (&c[i]);
            comp_inak (&c[i]);
            comp_ipca (&c[i]);
            comp_if (&c[i]);
            comp_istim (&c[i]);
            comp_itot (&c[i]);

            comp_ip3 (&c[i]);
            comp_qrel1 (&c[i]);
            comp_qrel2 (&c[i]);
            comp_qup1 (&c[i]);
            comp_qup2 (&c[i]);
            comp_qtr1 (&c[i]);
            comp_qtr2 (&c[i]);

            comp_conc (&c[i]);

            dvdt	    = -c[i].itot;
            c[i].v	   += dvdt*dt;

        }
        for (i = 0; i < n; i++) 
            printtofile(files[i],c[i]);
        t += dt;
    }
    free(files);
}

void printtofile (FILE *ap, Cell c)
{
	//count    += 1;
	fprintf(ap,"%f\t%f\t%f\n", t, c.v, c.caavg);
	//PRINT LAST 5 BEATS
	//if (count>=10 && t>=(BCL*(beats-5)))			
	//{
	//	count=0;
		//fprintf(ap,"%f\t%f\t%f\n", t-BCL*(beats-5), v, caavg);
   	//}	
	
}

void timestep (Cell *c, int n)
{
    int i;
    for (i = 0; i < n; i++)
        c[i].dvdtclock += dt;
}

void comp_conc (Cell *c)
{
	qdiff       = (c->casss-c->cassl)/sstau;  
	qgap        = (c->cassl-c->cai)/gaptau;  
    qdiffna     = (c->nasss-c->nassl)/sstau;
    qgapna      = (c->nassl-c->nai)/gaptau;

    //printf("qdiff = %lf\n",qdiff);
	//printf("qgap = %lf\n",qgap);
	//printf("qdiffna = %lf\n",qdiffna);
	//printf("qgapna = %lf\n",qgapna);
    
	dcasss		= dt*(-(c->ical-2*c->inacass)*acap/(vsss*2.0*frdy)+(c->qrel1+c->qip3)*vjsr/vsss-qdiff);
	bsss        = 1/(1+(bsrbar*kmbsr/pow(kmbsr+c->casss,2))+(bslbar*kmbsl/pow(kmbsl+c->casss,2)));
	c->casss      += bsss*dcasss;
    //printf("casss = %lf\n",c->casss);
	
	dcassl		= dt*(-(c->qup1)*vnsr/vssl+qdiff*vsss/vssl-qgap-(c->icat+c->ipca+c->icab-2*c->inaca)*acap/(vssl*2.0*frdy));
	trpn        = trpnbar1*(c->cassl/(c->cassl+kmtrpn));
	cmdn		= cmdnbar1*(c->cassl/(c->cassl+kmcmdn));
	catotal		= trpn+cmdn+dcassl+c->cassl;
	bmyo		= cmdnbar1+trpnbar1-catotal+kmtrpn+kmcmdn;
	cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar1*kmcmdn)+cmdnbar1*kmtrpn;
	dmyo		= -kmtrpn*kmcmdn*catotal;
	c->cassl		= (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;   
 	//printf("casss = %lf\n",c->cassl);

	dcajsr		= dt*(c->qtr1-c->qrel1-c->qip3);
	csqn1       = csqnbar1*(c->cajsr/(c->cajsr+kmcsqn));
	bjsr        = csqnbar1 - csqn1-c->cajsr-dcajsr+kmcsqn;
	cjsr        = kmcsqn*(csqn1+c->cajsr+dcajsr);
	c->cajsr       = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
	//printf("casss = %lf\n",c->cajsr);

	dcacsr		= dt*(c->qtr2-c->qrel2);
	csqn        = csqnbar*(c->cacsr/(c->cacsr+kmcsqn));
	bcsr        = csqnbar - csqn-c->cacsr-dcacsr+kmcsqn;
	ccsr        = kmcsqn*(csqn+c->cacsr+dcacsr);
	c->cacsr    = (sqrt(bcsr*bcsr+4*ccsr)-bcsr)/2;
	//printf("casss = %lf\n",c->cacsr);

	dcansr	    = dt*(c->qup1+c->qup2-c->qtr1*vjsr/vnsr-c->qtr2*vcsr/vnsr);
 	c->cansr	   += dcansr;
    //printf("cansr = %lf\n",c->cansr);
 	
	dnasss	    = dt*((-(3*c->inacass)*acap)/((vsss)*zna*frdy)-qdiffna); 
	c->nasss      += dnasss;
    //printf("nasss = %lf\n",c->nasss);
	
	dnassl	    = dt*((-(3*c->inak+c->ina+c->inal+3*c->inaca+c->ifna+c->inab)*acap)/((vssl)*zna*frdy)+qdiffna*vsss/vssl-qgapna);
	c->nassl	   += dnassl;
    //printf("nassl = %lf\n",c->nassl);
	
	dnai        = dt*(qgapna*vssl/vmyo);
	c->nai        += dnai;
    //printf("nai = %lf\n",c->nai);
	
	dki	        = dt*((-c->iktot*acap)/((vmyo+vssl+vsss)*zk*frdy));
	c->ki         += dki;
    //printf("ki = %lf\n",c->ki);
	
	dcai		= dt*(-(c->qup2)*vnsr/vmyo+qgap*vssl/vmyo+(c->qrel2)*vcsr/vmyo);
	trpn        = trpnbar*(c->cai/(c->cai+kmtrpn));
	cmdn		= cmdnbar*(c->cai/(c->cai+kmcmdn));
	catotal		= trpn+cmdn+dcai+c->cai;
	bmyo		= cmdnbar+trpnbar-catotal+kmtrpn+kmcmdn;
	cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar*kmcmdn)+cmdnbar*kmtrpn;
	dmyo		= -kmtrpn*kmcmdn*catotal;
	c->cai		    = (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;  
	
	c->caavg       = (c->casss*vsss+c->cassl*vssl+c->cai*vmyo)/(vsss+vmyo+vssl);
	//printf("caavg = %lf\n",c->caavg);

 	camkbound	= camk0*(1-c->camktrap)*1/(1+(kmcam/c->casss));
	c->camktrap	= dt*(alphacamk*camkbound*(camkbound+c->camktrap)-betacamk*c->camktrap) + c->camktrap;
	c->camkactive	= camkbound+c->camktrap; 
	//printf("camkactive = %lf\n",c->camkactive);

}         

void comp_qtr2 (Cell *c)
{
	c->qtr2		= (c->cansr-c->cacsr)/tautr2;
}

void comp_qtr1 (Cell *c)
{
	c->qtr1		= (c->cansr-c->cajsr)/tautr1;
}

void comp_qup2 (Cell *c)
{
    dkmplb		= dkmplbbar*c->camkactive/(kmcamk+c->camkactive);
	dqupcamk	= dqupcamkbar*c->camkactive/(kmcamk+c->camkactive); 
	c->qup2		= 0.0026*(dqupcamk+1)/(1+pow((kmup-dkmplb)/c->cai,1))-0.0042*c->cansr/nsrbar;
}

void comp_qup1 (Cell *c)
{
    dkmplb		= dkmplbbar*c->camkactive/(kmcamk+c->camkactive);
	dqupcamk	= dqupcamkbar*c->camkactive/(kmcamk+c->camkactive); 
	c->qup1		= 0.0002*(dqupcamk+1)/(1+pow((kmup-dkmplb)/c->cassl,1))-0.00105*c->cansr/nsrbar;
}

void comp_qrel2 (Cell *c)
{
	qgap  = (c->cassl-c->cai)/gaptau;  
    REL  = (-c->qup2*vnsr/vmyo + qgap*vssl/vmyo+ (c->qrel2)*vcsr/vmyo);    
    ireltau = 6*(1+1*(1/(1+pow((0.28/c->camkactive),8))))/(1+(0.0123/c->cacsr));
    if (REL > 0)
        irelss  = 91*(1+1*(1/(1+pow((0.28/c->camkactive),8))))*(REL)/(1 + pow((1/c->cacsr),8));
    else 
        irelss = 0;
    c->qrel2 += dt*((irelss-c->qrel2)/ireltau);
}

void comp_qrel1 (Cell *c)
{
	qdiff  = (c->casss-c->cassl)/sstau;  
    REL  = -((c->ical)*acap/(vsss*2.0*frdy) - (c->qrel1 + c->qip3)*vjsr/vsss + qdiff);     
    ireltau = 2*(1+1*(1/(1+pow((0.28/c->camkactive),8))))/(1+(0.0123/c->cajsr));
    if (REL > 0)
        irelss  = 15*(1+1*(1/(1+pow((0.28/c->camkactive),8))))*REL/(1 + pow((1.0/c->cajsr),8));
    else 
        irelss = 0;
    c->qrel1 += dt*((irelss-c->qrel1)/ireltau);
}

void comp_ip3 (Cell *c)
{
    c->u += dt*(c->casss*k2*(1-c->u) - k2a*c->u);
    POip3 = tauip3r*IP3*c->casss*(1-c->u)/((1+IP3*k0/k0a)*(1+c->casss*k1/k1a));
    c->qip3 = 10.920*(c->cajsr-c->casss)*(POip3);
}

void comp_itot (Cell *c)
{
	if (c->stimtime >= 0.0 && c->stimtime < stimdur)
	{
		//printf("stimtime = %lf || stimdur = %lf\n",stimtime,stimdur);
		c->icatot	= c->ical+c->icat+c->ipca+c->icab-2*c->inaca-2*c->inacass;
		c->iktot	= c->ikr+c->iks+c->ik1-2*c->inak+c->ito1+c->ifk+1*istim;
		c->inatot	= 3*c->inak+c->ina+3*c->inaca+3*c->inacass+c->inal+c->ifna+c->inab;
		c->itot	= c->icatot+c->iktot+c->inatot;
	}
	else
	{
		c->icatot	= c->ical+c->icat+c->ipca+c->icab-2*c->inaca-2*c->inacass;
		c->iktot	= c->ikr+c->iks+c->ik1-2*c->inak+c->ito1+c->ifk;
		c->inatot	= 3*c->inak+c->ina+3*c->inaca+3*c->inacass+c->inal+c->ifna+c->inab;
		c->itot	= c->icatot+c->iktot+c->inatot;
	}
}

void comp_istim (Cell *c) 
{
	c->stimtime += dt;
	if (t >= c->tstim)
	{
		c->stimtime = 0.0;					
		c->stimcount += 1;					
		if (c->stimcount < c->beats-1)  c->tstim += c->BCL;		
		else if (c->stimcount == c->beats-1) c->tstim += c->S2;	
		else c->tstim = tmax+1;				
		//if (c->stimcount < c->beats) printf ("S1 Beat %d at time = %.2f ms !\n", stimcount+1, t);
		//else if (c->stimcount == c->beats) printf ("S2 Beat at time = %.2f ms !\n", t);
	}
}

void comp_if (Cell *c)
{
	yss       = 1/(1+exp((c->v+87)/9.5));
	ytau      = 2000/(exp(-(c->v+132)/10) + exp((c->v+57)/60));
	c->y      = yss - (yss-c->y)*exp(-dt/ytau);
	c->ifna	  = 0.012*c->y*c->y*(c->v-ena);
	c->ifk       = 0.024*c->y*c->y*(c->v-ek);
	c->iftotal   = c->ifna + c->ifk;
}

void comp_ipca (Cell *c)
{
    c->ipca = ipcabar/((kmpca/c->cassl)+1);
}

void comp_inak (Cell *c)
{
    c->inak	= ibarnak*(1/(1+exp(-1*(c->v+92)*frdy/(R*temp))))*pow((c->nassl/(c->nassl+2.6)),3)*(ko/(ko+0.8));
}

void comp_inaca (Cell *c)
{
	allo		= 1/(1+pow((kmcaact/(1.5*c->casss)),2));
	num		    = inacamax*(pow(c->nasss,3)*cao*exp(nu*c->v*frdy/(R*temp))-pow(nao,3)*1.5*c->casss*exp((nu-1)*c->v*frdy/(R*temp)));
	denommult	= 1+ksat*exp((nu-1)*c->v*frdy/(R*temp));
	denomterm1	= kmcao*pow(c->nasss,3)+pow(kmnao,3)*1.5*c->casss+pow(kmnai1,3)*cao*(1+1.5*c->casss/kmcai);
	denomterm2	= kmcai*pow(nao,3)*(1+pow(c->nasss/kmnai1,3))+pow(c->nasss,3)*cao+pow(nao,3)*1.5*c->casss;
	deltaE		= num/(denommult*(denomterm1+denomterm2));
	c->inacass  = 0.2*allo*deltaE;
	
	allo		= 1/(1+pow((kmcaact/(1.5*c->cassl)),2));
	num		    = inacamax*(pow(c->nassl,3)*cao*exp(nu*c->v*frdy/(R*temp))-pow(nao,3)*1.5*c->cassl*exp((nu-1)*c->v*frdy/(R*temp)));
	denommult	= 1+ksat*exp((nu-1)*c->v*frdy/(R*temp));
	denomterm1	= kmcao*pow(c->nassl,3)+pow(kmnao,3)*1.5*c->cassl+pow(kmnai1,3)*cao*(1+1.5*c->cassl/kmcai);
	denomterm2	= kmcai*pow(nao,3)*(1+pow(c->nassl/kmnai1,3))+pow(c->nassl,3)*cao+pow(nao,3)*1.5*c->cassl;
	deltaE		= num/(denommult*(denomterm1+denomterm2));
	c->inaca    = 0.8*allo*deltaE;
}

void comp_ik1 (Cell *c)
{	
    k1ss      = 1/(1+exp((c->v+103-(2.9+ko*2.175))/10.15));
	gk1	      = 0.12*sqrt(ko);
	c->ik1	      = gk1*k1ss*(c->v-ek);
}

void comp_iks (Cell *c)
{
	eks	    = (R*temp/frdy)*log((ko+prnak*nao)/(c->ki+prnak*c->nassl));
	gks	    = 0.053*(1+0.6/(1+pow((0.000038/c->cassl),1.4)));
	xsss	= 1/(1+exp(-(c->v-9)/13.7));
	xs1tau	= 200/(exp(-(c->v+10)/6) + exp((c->v-62)/55));
	xs2tau	= 1500+ 350/(exp(-(c->v+10)/4) + exp((c->v-90)/58));
	c->xs1	= xsss-(xsss-c->xs1)*exp(-dt/xs1tau);
	c->xs2	= xsss-(xsss-c->xs2)*exp(-dt/xs2tau);
	c->iks	= gks*c->xs1*c->xs2*(c->v-eks);
}

void comp_ikr (Cell *c)
{
	gkr	    = 0.0326*sqrt(ko/5.4);
	xrss	= 1/(1+exp(-(c->v)/15));
	xrtau   = 400.0/(1.0+exp(c->v/10.0)) + 100.0;
	rkr	    = 1/(1+exp((c->v)/35));
	c->xr	    = xrss-(xrss-c->xr)*exp(-dt/xrtau);
	c->ikr	    = gkr*c->xr*rkr*(c->v-ek);
}

void comp_ito1 (Cell *c)
{
	atau	= 1/(25*exp((c->v-82)/18)/(1+exp((c->v-82)/18))+25*exp(-(c->v+52)/18)/(1+exp(-(c->v+52)/18)));
	itau	= 2.86+ 1/(exp(-(c->v+125)/15)*0.1 + 0.1*exp((c->v+2)/26.5));
	i2tau	= 21.5+ 1/(exp(-(c->v+138.2)/52)*0.005 + 0.003*exp((c->v+18)/12.5));
	ass	    = 1/(1+exp(-(c->v-8.9)/10.3));
	iss	    = 1/(1+exp((c->v+30)/11));
	i2ss	= iss;
	c->a	    = ass-(ass-c->a)*exp(-dt/atau);
	c->i	    = iss-(iss-c->i)*exp(-dt/itau);
	c->i2	    = i2ss-(i2ss-c->i2)*exp(-dt/i2tau);
	c->itos    = gtos*c->a*c->i*c->i2*(c->v-ek);
	c->itof    = gtof*(c->v-ek)/(1+exp(-(c->v-3)/19.8));
	c->ito1	= c->itos + c->itof;
}

void comp_icab (Cell *c)
{
	c->icab	= pcab*zca*zca*((c->v*frdy*frdy)/(R*temp))*((gacai*c->cassl*exp((zca*c->v*frdy)/(R*temp))-gacao*cao)/(exp((zca*c->v*frdy)/(R*temp))-1));
}

void comp_icat (Cell *c)
{
	bss	    = 1/(1+ exp (-(c->v+30)/7));
	gss	    = 1/(1+exp((c->v+61)/5));
	taub	= 1/(1.068*exp((c->v+16.3)/30)+1.068*exp(-(c->v+16.3)/30));
	taug    = 1/(0.015*exp(-(c->v+71.7)/83.3)+0.015*exp((c->v+71.7)/15.4));
	c->b	= bss-(bss-c->b)*exp(-dt/taub);
	c->g	= gss-(gss-c->g)*exp(-dt/taug);
	c->icat	= gcat*c->b*c->g*(c->v-eca);
}

/* Pode ser que precise colocar as corrente 'ical' como variavel da struct Cell */
void comp_ical (Cell *c)
{
	ibarca		= pca*zca*zca*(((c->v-15)*frdy*frdy)/(R*temp))*((gacai*c->casss*exp((zca*(c->v-15)*frdy)/(R*temp))-gacao*cao)/(exp((zca*(c->v-15)*frdy)/(R*temp))-1));
	dss		    = (1/(1.0+exp(-(c->v-2.0)/7.8)));
	dtau		= (0.59+0.8*exp(0.052*(c->v+13))/(1+exp(0.132*(c->v+13))));
	fss	        = 1/(1.0 + exp((c->v+16.5)/9.5));
	ftau        = 0.92/(0.125*exp(-(0.058*(c->v-2.5))*(0.045*(c->v-2.5)))+0.1);
	f2ss        = fss;
	f2tau       = 0.90/(0.02*exp(-(0.04*(c->v-18.6))*(0.045*(c->v-18.6)))+0.005);
	fcass		= 0.3/(1 - c->ical/0.05) + 0.55/(1.0+c->casss/0.003)+0.15;
	fcatau		= 10*c->camkactive/(c->camkactive+kmcam) + 0.5+1/(1.0+c->casss/0.003);
	fca2ss		= 1.0/(1.0-c->ical/0.01);
	fca2tau		= 1*(300.0/(1.0+exp((-c->ical-0.175)/0.04))+125.0);
	c->d		= dss-(dss-c->d)*exp(-dt/dtau);
	c->f		= fss-(fss-c->f)*exp(-dt/ftau);
	c->f2		= f2ss-(f2ss-c->f2)*exp(-dt/f2tau);
	c->fca		= fcass-(fcass-c->fca)*exp(-dt/fcatau);
	c->fca2		= fca2ss-(fca2ss-c->fca2)*exp(-dt/fca2tau);
	c->ical		= c->d*c->f*c->f2*c->fca*c->fca2*ibarca;	
}

void comp_inab (Cell *c)
{
    c->inab    = pnab*frdy*((frdy*c->v)/(R*temp))*(c->nassl*exp((frdy*c->v)/(R*temp)) - nao)/(exp((frdy*c->v)/(R*temp))-1);     
}

void comp_inal (Cell *c)
{
	mltau	= 1/(0.64*(c->v+37.13)/(1-exp(-0.1*(c->v+37.13))) + 0.16*exp(-c->v/11));
	ml3tau  = mltau;
	mlss	= 1/(1+exp(-(c->v+28)/7));
	ml3ss   = 1/(1+exp(-(c->v+63)/7));
	hltau   = 162+132/(1+exp(-(c->v+28)/5.5));
	hl3tau  = 0.5*hltau;
	hlss	= 1/(1+exp((c->v+28)/12));
	hl3ss	= 1/(1+exp((c->v+63)/12));
	jltau   = 411;
	jl3tau  = 0.5*jltau;
	jlss	= hlss;
	jl3ss	= hl3ss;
	c->ml	  = mlss-(mlss-c->ml)*exp(-dt/mltau);
	c->ml3     = ml3ss-(ml3ss-c->ml3)*exp(-dt/ml3tau);
	c->hl	  = hlss-(hlss-c->hl)*exp(-dt/hltau);
	c->hl3     = hl3ss-(hl3ss-c->hl3)*exp(-dt/hl3tau);
	c->jl	  = jlss-(jlss-c->jl)*exp(-dt/jltau);
	c->jl3     = jl3ss-(jl3ss-c->jl3)*exp(-dt/jl3tau);
	c->inal2   = gnal2*c->ml*c->hl*c->jl*(c->v-ena);
	c->inal3   = gnal3*c->ml3*c->hl3*c->jl3*(c->v-ena);
	c->inal    = c->inal2 + c->inal3; 
}

void comp_ina (Cell *c)
{
    ma	= 0.64*(c->v+37.13)/(1-exp(-0.1*(c->v+37.13)));
	mb	= 0.16*exp(-c->v/11);
	if (c->v<-40)
	{
		ha = 0.135*exp((70+c->v)/-6.8);
		hb = 3.56*exp(0.079*c->v)+310000*exp(0.35*c->v);
		ja = (-127140*exp(0.2444*c->v)-0.003474*exp(-0.04391*c->v))*(c->v+37.78)/(1+exp(0.311*(c->v+79.23)));
		jb = 0.1212*exp(-0.01052*c->v)/(1+exp(-0.1378*(c->v+40.14)));
	}
	else
	{
		ha = 0.0;
		hb = 1/(0.13*(1+exp((c->v+10.66)/-11.1)));
		ja = 0.0;
		jb = 0.3*exp(-0.0000002535*c->v)/(1+exp(-0.1*(c->v+32)));
	}
	mtau	= 1/(ma+mb);
	htau	= 1/(ha+hb);
	jtau	= 1/(ja+jb);
	mss	= ma*mtau;
	hss	= ha*htau;
	jss	= 1*ja*jtau;
	c->m	= mss-(mss-c->m)*exp(-dt/mtau);
	c->h	= hss-(hss-c->h)*exp(-dt/htau);
	c->j	= jss-(jss-c->j)*exp(-dt/jtau);
	c->ina	= gna*pow(c->m,3)*c->h*c->j*(c->v-ena);
}

void comp_revs (Cell *c)
{
    eca	= (R*temp/(zca*frdy))*log(cao/c->cassl);
	ena	= (R*temp/frdy)*log(nao/c->nassl);
	ek	= (R*temp/frdy)*log(ko/c->ki);
}

FILE** createFiles (int n)
{
    int i;
    FILE **files = (FILE**)malloc(sizeof(FILE*)*n);
    char filename[MAX_FILENAME];
    for (i = 0; i < n; i++)
    {   
        sprintf(filename,"cell%d.dat",i);
        files[i] = fopen(filename,"w+");
    }
    return files;
}

void printCells (Cell c[], int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        printf("---- Cell %d -----\n",i);
        printf("Beats = %d\n",c[i].beats);
        printf("BCL = %lf\n",c[i].BCL);
        printf("S2 = %lf\n",c[i].S2);
        printf("tstim = %lf\n",c[i].tstim);
    }
}