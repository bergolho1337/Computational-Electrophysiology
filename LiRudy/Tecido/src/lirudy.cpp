#include "../include/lirudy.h"

// TIME
double t, T;

// VOLTAGE
double v,dvdtclock;

// STIMULUS CURRENT
double istim = -80;	        
double tstim,stimtime;			   
int stimcount;

// CALCIUM FLUXES AND CONCENTRATIONS	
double IP3 = 0.0001;

double kmup   = 0.00028;

double trpnbar1 = 3.5e-3;
double cmdnbar1 = 1.25e-2;
double csqnbar1 = 1.2;

// OUTPUT FILE
int count;

LiRudy::LiRudy (int argc, char *argv[])
{
	dt = atof(argv[1]);
	tmax = atof(argv[2]);
	n = tmax / dt;

	computeGeometrics();
	allocMem();
	setInitCond();
	setStimCells();
	


}

void LiRudy::setInitCond ()
{
	for (int i = 0; i < (int)cells.size(); i++)
    	cells[i].setCell();
}

void LiRudy::setStimCells ()
{
	for (int i = 0; i < NSC; i++)
    {
        cells[i].beats = 10;
        cells[i].BCL = 500;
        cells[i].S2 = 500;
        cells[i].tstim = 0;
    }
    for (int i = NSC; i < NCELLS; i++)
    {
        cells[i].beats = 10;
        cells[i].BCL = tmax + 500;
        cells[i].S2 = 500;
        cells[i].tstim = cells[i].BCL;
    }
}

void LiRudy::allocMem ()
{
	cells.assign(NCELLS,Cell());
}

/*

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
	//#pragma omp parallel for num_threads(2)
	for (int i = 0; i < n; i++)
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
    FILE **files = createFiles(1);
    int i;
    // Time loop
    while (t <= tmax)
    {
        timestep(c,n);
        // Space loop
		//#pragma omp parallel for num_threads(2)
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

			c[i].dvdt	   = -c[i].itot;
            c[i].v	   += c[i].dvdt*dt;

        }
        //for (i = 0; i < n; i++) 
            printtofile(files[0],c[0]);
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
	c->qdiff       = (c->casss-c->cassl)/sstau;  
	c->qgap        = (c->cassl-c->cai)/gaptau;  
    c->qdiffna     = (c->nasss-c->nassl)/sstau;
    c->qgapna      = (c->nassl-c->nai)/gaptau;

    //printf("qdiff = %lf\n",qdiff);
	//printf("qgap = %lf\n",qgap);
	//printf("qdiffna = %lf\n",qdiffna);
	//printf("qgapna = %lf\n",qgapna);
    
	c->dcasss		= dt*(-(c->ical-2*c->inacass)*acap/(vsss*2.0*frdy)+(c->qrel1+c->qip3)*vjsr/vsss-c->qdiff);
	c->bsss        = 1/(1+(bsrbar*kmbsr/pow(kmbsr+c->casss,2))+(bslbar*kmbsl/pow(kmbsl+c->casss,2)));
	c->casss      += c->bsss*c->dcasss;
    //printf("casss = %lf\n",c->casss);
	
	c->dcassl		= dt*(-(c->qup1)*vnsr/vssl+c->qdiff*vsss/vssl-c->qgap-(c->icat+c->ipca+c->icab-2*c->inaca)*acap/(vssl*2.0*frdy));
	c->trpn        = trpnbar1*(c->cassl/(c->cassl+kmtrpn));
	c->cmdn		= cmdnbar1*(c->cassl/(c->cassl+kmcmdn));
	c->catotal		= c->trpn+c->cmdn+c->dcassl+c->cassl;
	c->bmyo		= cmdnbar1+trpnbar1-c->catotal+kmtrpn+kmcmdn;
	c->cmyo		= kmcmdn*kmtrpn-c->catotal*(kmtrpn+kmcmdn)+(trpnbar1*kmcmdn)+cmdnbar1*kmtrpn;
	c->dmyo		= -kmtrpn*kmcmdn*c->catotal;
	c->cassl		= (2.0/3.0)*sqrt(c->bmyo*c->bmyo-3.0*c->cmyo)*cos(acos((9.0*c->bmyo*c->cmyo-2*c->bmyo*c->bmyo*c->bmyo-27*c->dmyo)/(2.0*pow((c->bmyo*c->bmyo-3.0*c->cmyo),1.5)))/3.0)-c->bmyo/3.0;   
 	//printf("casss = %lf\n",c->cassl);

	c->dcajsr		= dt*(c->qtr1-c->qrel1-c->qip3);
	c->csqn1       = csqnbar1*(c->cajsr/(c->cajsr+kmcsqn));
	c->bjsr        = csqnbar1 - c->csqn1-c->cajsr-c->dcajsr+kmcsqn;
	c->cjsr        = kmcsqn*(c->csqn1+c->cajsr+c->dcajsr);
	c->cajsr       = (sqrt(c->bjsr*c->bjsr+4*c->cjsr)-c->bjsr)/2;
	//printf("casss = %lf\n",c->cajsr);

	c->dcacsr		= dt*(c->qtr2-c->qrel2);
	c->csqn        = csqnbar*(c->cacsr/(c->cacsr+kmcsqn));
	c->bcsr        = csqnbar - c->csqn-c->cacsr-c->dcacsr+kmcsqn;
	c->ccsr        = kmcsqn*(c->csqn+c->cacsr+c->dcacsr);
	c->cacsr    = (sqrt(c->bcsr*c->bcsr+4*c->ccsr)-c->bcsr)/2;
	//printf("casss = %lf\n",c->cacsr);

	c->dcansr	    = dt*(c->qup1+c->qup2-c->qtr1*vjsr/vnsr-c->qtr2*vcsr/vnsr);
 	c->cansr	   += c->dcansr;
    //printf("cansr = %lf\n",c->cansr);
 	
	c->dnasss	    = dt*((-(3*c->inacass)*acap)/((vsss)*zna*frdy)-c->qdiffna); 
	c->nasss      += c->dnasss;
    //printf("nasss = %lf\n",c->nasss);
	
	c->dnassl	    = dt*((-(3*c->inak+c->ina+c->inal+3*c->inaca+c->ifna+c->inab)*acap)/((vssl)*zna*frdy)+c->qdiffna*vsss/vssl-c->qgapna);
	c->nassl	   += c->dnassl;
    //printf("nassl = %lf\n",c->nassl);
	
	c->dnai        = dt*(c->qgapna*vssl/vmyo);
	c->nai        += c->dnai;
    //printf("nai = %lf\n",c->nai);
	
	c->dki	        = dt*((-c->iktot*acap)/((vmyo+vssl+vsss)*zk*frdy));
	c->ki         += c->dki;
    //printf("ki = %lf\n",c->ki);
	
	c->dcai		= dt*(-(c->qup2)*vnsr/vmyo+c->qgap*vssl/vmyo+(c->qrel2)*vcsr/vmyo);
	c->trpn        = trpnbar*(c->cai/(c->cai+kmtrpn));
	c->cmdn		= cmdnbar*(c->cai/(c->cai+kmcmdn));
	c->catotal		= c->trpn+c->cmdn+c->dcai+c->cai;
	c->bmyo		= cmdnbar+trpnbar-c->catotal+kmtrpn+kmcmdn;
	c->cmyo		= kmcmdn*kmtrpn-c->catotal*(kmtrpn+kmcmdn)+(trpnbar*kmcmdn)+cmdnbar*kmtrpn;
	c->dmyo		= -kmtrpn*kmcmdn*c->catotal;
	c->cai		    = (2.0/3.0)*sqrt(c->bmyo*c->bmyo-3.0*c->cmyo)*cos(acos((9.0*c->bmyo*c->cmyo-2*c->bmyo*c->bmyo*c->bmyo-27*c->dmyo)/(2.0*pow((c->bmyo*c->bmyo-3.0*c->cmyo),1.5)))/3.0)-c->bmyo/3.0;  
	
	c->caavg       = (c->casss*vsss+c->cassl*vssl+c->cai*vmyo)/(vsss+vmyo+vssl);
	//printf("caavg = %lf\n",c->caavg);

 	c->camkbound	= camk0*(1-c->camktrap)*1/(1+(kmcam/c->casss));
	c->camktrap	= dt*(alphacamk*c->camkbound*(c->camkbound+c->camktrap)-betacamk*c->camktrap) + c->camktrap;
	c->camkactive	= c->camkbound+c->camktrap; 
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
    c->dkmplb		= dkmplbbar*c->camkactive/(kmcamk+c->camkactive);
	c->dqupcamk	= dqupcamkbar*c->camkactive/(kmcamk+c->camkactive); 
	c->qup2		= 0.0026*(c->dqupcamk+1)/(1+pow((kmup-c->dkmplb)/c->cai,1))-0.0042*c->cansr/nsrbar;
}

void comp_qup1 (Cell *c)
{
    c->dkmplb		= dkmplbbar*c->camkactive/(kmcamk+c->camkactive);
	c->dqupcamk	= dqupcamkbar*c->camkactive/(kmcamk+c->camkactive); 
	c->qup1		= 0.0002*(c->dqupcamk+1)/(1+pow((kmup-c->dkmplb)/c->cassl,1))-0.00105*c->cansr/nsrbar;
}

void comp_qrel2 (Cell *c)
{
	c->qgap  = (c->cassl-c->cai)/gaptau;  
    c->REL  = (-c->qup2*vnsr/vmyo + c->qgap*vssl/vmyo+ (c->qrel2)*vcsr/vmyo);    
    c->ireltau = 6*(1+1*(1/(1+pow((0.28/c->camkactive),8))))/(1+(0.0123/c->cacsr));
    if (c->REL > 0)
        c->irelss  = 91*(1+1*(1/(1+pow((0.28/c->camkactive),8))))*(c->REL)/(1 + pow((1/c->cacsr),8));
    else 
        c->irelss = 0;
    c->qrel2 += dt*((c->irelss-c->qrel2)/c->ireltau);
}

void comp_qrel1 (Cell *c)
{
	c->qdiff  = (c->casss-c->cassl)/sstau;  
    c->REL  = -((c->ical)*acap/(vsss*2.0*frdy) - (c->qrel1 + c->qip3)*vjsr/vsss + c->qdiff);     
    c->ireltau = 2*(1+1*(1/(1+pow((0.28/c->camkactive),8))))/(1+(0.0123/c->cajsr));
    if (c->REL > 0)
        c->irelss  = 15*(1+1*(1/(1+pow((0.28/c->camkactive),8))))*c->REL/(1 + pow((1.0/c->cajsr),8));
    else 
        c->irelss = 0;
    c->qrel1 += dt*((c->irelss-c->qrel1)/c->ireltau);
}

void comp_ip3 (Cell *c)
{
    c->u += dt*(c->casss*k2*(1-c->u) - k2a*c->u);
    c->POip3 = tauip3r*IP3*c->casss*(1-c->u)/((1+IP3*k0/k0a)*(1+c->casss*k1/k1a));
    c->qip3 = 10.920*(c->cajsr-c->casss)*(c->POip3);
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
	c->yss       = 1/(1+exp((c->v+87)/9.5));
	c->ytau      = 2000/(exp(-(c->v+132)/10) + exp((c->v+57)/60));
	c->y      = c->yss - (c->yss-c->y)*exp(-dt/c->ytau);
	c->ifna	  = 0.012*c->y*c->y*(c->v-c->ena);
	c->ifk       = 0.024*c->y*c->y*(c->v-c->ek);
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
	c->allo		= 1/(1+pow((kmcaact/(1.5*c->casss)),2));
	c->num		    = inacamax*(pow(c->nasss,3)*cao*exp(nu*c->v*frdy/(R*temp))-pow(nao,3)*1.5*c->casss*exp((nu-1)*c->v*frdy/(R*temp)));
	c->denommult	= 1+ksat*exp((nu-1)*c->v*frdy/(R*temp));
	c->denomterm1	= kmcao*pow(c->nasss,3)+pow(kmnao,3)*1.5*c->casss+pow(kmnai1,3)*cao*(1+1.5*c->casss/kmcai);
	c->denomterm2	= kmcai*pow(nao,3)*(1+pow(c->nasss/kmnai1,3))+pow(c->nasss,3)*cao+pow(nao,3)*1.5*c->casss;
	c->deltaE		= c->num/(c->denommult*(c->denomterm1+c->denomterm2));
	c->inacass  = 0.2*c->allo*c->deltaE;
	
	c->allo		= 1/(1+pow((kmcaact/(1.5*c->cassl)),2));
	c->num		    = inacamax*(pow(c->nassl,3)*cao*exp(nu*c->v*frdy/(R*temp))-pow(nao,3)*1.5*c->cassl*exp((nu-1)*c->v*frdy/(R*temp)));
	c->denommult	= 1+ksat*exp((nu-1)*c->v*frdy/(R*temp));
	c->denomterm1	= kmcao*pow(c->nassl,3)+pow(kmnao,3)*1.5*c->cassl+pow(kmnai1,3)*cao*(1+1.5*c->cassl/kmcai);
	c->denomterm2	= kmcai*pow(nao,3)*(1+pow(c->nassl/kmnai1,3))+pow(c->nassl,3)*cao+pow(nao,3)*1.5*c->cassl;
	c->deltaE		= c->num/(c->denommult*(c->denomterm1+c->denomterm2));
	c->inaca    = 0.8*c->allo*c->deltaE;
}

void comp_ik1 (Cell *c)
{	
    c->k1ss      = 1/(1+exp((c->v+103-(2.9+ko*2.175))/10.15));
	c->gk1	      = 0.12*sqrt(ko);
	c->ik1	      = c->gk1*c->k1ss*(c->v-c->ek);
}

void comp_iks (Cell *c)
{
	c->eks	    = (R*temp/frdy)*log((ko+prnak*nao)/(c->ki+prnak*c->nassl));
	c->gks	    = 0.053*(1+0.6/(1+pow((0.000038/c->cassl),1.4)));
	c->xsss	= 1/(1+exp(-(c->v-9)/13.7));
	c->xs1tau	= 200/(exp(-(c->v+10)/6) + exp((c->v-62)/55));
	c->xs2tau	= 1500+ 350/(exp(-(c->v+10)/4) + exp((c->v-90)/58));
	c->xs1	= c->xsss-(c->xsss-c->xs1)*exp(-dt/c->xs1tau);
	c->xs2	= c->xsss-(c->xsss-c->xs2)*exp(-dt/c->xs2tau);
	c->iks	= c->gks*c->xs1*c->xs2*(c->v-c->eks);
}

void comp_ikr (Cell *c)
{
	c->gkr	    = 0.0326*sqrt(ko/5.4);
	c->xrss	= 1/(1+exp(-(c->v)/15));
	c->xrtau   = 400.0/(1.0+exp(c->v/10.0)) + 100.0;
	c->rkr	    = 1/(1+exp((c->v)/35));
	c->xr	    = c->xrss-(c->xrss-c->xr)*exp(-dt/c->xrtau);
	c->ikr	    = c->gkr*c->xr*c->rkr*(c->v-c->ek);
}

void comp_ito1 (Cell *c)
{
	c->atau	= 1/(25*exp((c->v-82)/18)/(1+exp((c->v-82)/18))+25*exp(-(c->v+52)/18)/(1+exp(-(c->v+52)/18)));
	c->itau	= 2.86+ 1/(exp(-(c->v+125)/15)*0.1 + 0.1*exp((c->v+2)/26.5));
	c->i2tau	= 21.5+ 1/(exp(-(c->v+138.2)/52)*0.005 + 0.003*exp((c->v+18)/12.5));
	c->ass	    = 1/(1+exp(-(c->v-8.9)/10.3));
	c->iss	    = 1/(1+exp((c->v+30)/11));
	c->i2ss	= c->iss;
	c->a	    = c->ass-(c->ass-c->a)*exp(-dt/c->atau);
	c->i	    = c->iss-(c->iss-c->i)*exp(-dt/c->itau);
	c->i2	    = c->i2ss-(c->i2ss-c->i2)*exp(-dt/c->i2tau);
	c->itos    = gtos*c->a*c->i*c->i2*(c->v-c->ek);
	c->itof    = gtof*(c->v-c->ek)/(1+exp(-(c->v-3)/19.8));
	c->ito1	= c->itos + c->itof;
}

void comp_icab (Cell *c)
{
	c->icab	= pcab*zca*zca*((c->v*frdy*frdy)/(R*temp))*((gacai*c->cassl*exp((zca*c->v*frdy)/(R*temp))-gacao*cao)/(exp((zca*c->v*frdy)/(R*temp))-1));
}

void comp_icat (Cell *c)
{
	c->bss	    = 1/(1+ exp (-(c->v+30)/7));
	c->gss	    = 1/(1+exp((c->v+61)/5));
	c->taub	= 1/(1.068*exp((c->v+16.3)/30)+1.068*exp(-(c->v+16.3)/30));
	c->taug    = 1/(0.015*exp(-(c->v+71.7)/83.3)+0.015*exp((c->v+71.7)/15.4));
	c->b	= c->bss-(c->bss-c->b)*exp(-dt/c->taub);
	c->g	= c->gss-(c->gss-c->g)*exp(-dt/c->taug);
	c->icat	= gcat*c->b*c->g*(c->v-c->eca);
}

// Pode ser que precise colocar as corrente 'ical' como variavel da struct Cell 
void comp_ical (Cell *c)
{
	c->ibarca		= pca*zca*zca*(((c->v-15)*frdy*frdy)/(R*temp))*((gacai*c->casss*exp((zca*(c->v-15)*frdy)/(R*temp))-gacao*cao)/(exp((zca*(c->v-15)*frdy)/(R*temp))-1));
	c->dss		    = (1/(1.0+exp(-(c->v-2.0)/7.8)));
	c->dtau		= (0.59+0.8*exp(0.052*(c->v+13))/(1+exp(0.132*(c->v+13))));
	c->fss	        = 1/(1.0 + exp((c->v+16.5)/9.5));
	c->ftau        = 0.92/(0.125*exp(-(0.058*(c->v-2.5))*(0.045*(c->v-2.5)))+0.1);
	c->f2ss        = c->fss;
	c->f2tau       = 0.90/(0.02*exp(-(0.04*(c->v-18.6))*(0.045*(c->v-18.6)))+0.005);
	c->fcass		= 0.3/(1 - c->ical/0.05) + 0.55/(1.0+c->casss/0.003)+0.15;
	c->fcatau		= 10*c->camkactive/(c->camkactive+kmcam) + 0.5+1/(1.0+c->casss/0.003);
	c->fca2ss		= 1.0/(1.0-c->ical/0.01);
	c->fca2tau		= 1*(300.0/(1.0+exp((-c->ical-0.175)/0.04))+125.0);
	c->d		= c->dss-(c->dss-c->d)*exp(-dt/c->dtau);
	c->f		= c->fss-(c->fss-c->f)*exp(-dt/c->ftau);
	c->f2		= c->f2ss-(c->f2ss-c->f2)*exp(-dt/c->f2tau);
	c->fca		= c->fcass-(c->fcass-c->fca)*exp(-dt/c->fcatau);
	c->fca2		= c->fca2ss-(c->fca2ss-c->fca2)*exp(-dt/c->fca2tau);
	c->ical		= c->d*c->f*c->f2*c->fca*c->fca2*c->ibarca;	
}

void comp_inab (Cell *c)
{
    c->inab    = pnab*frdy*((frdy*c->v)/(R*temp))*(c->nassl*exp((frdy*c->v)/(R*temp)) - nao)/(exp((frdy*c->v)/(R*temp))-1);     
}

void comp_inal (Cell *c)
{
	c->mltau	= 1/(0.64*(c->v+37.13)/(1-exp(-0.1*(c->v+37.13))) + 0.16*exp(-c->v/11));
	c->ml3tau  = c->mltau;
	c->mlss	= 1/(1+exp(-(c->v+28)/7));
	c->ml3ss   = 1/(1+exp(-(c->v+63)/7));
	c->hltau   = 162+132/(1+exp(-(c->v+28)/5.5));
	c->hl3tau  = 0.5*c->hltau;
	c->hlss	= 1/(1+exp((c->v+28)/12));
	c->hl3ss	= 1/(1+exp((c->v+63)/12));
	c->jltau   = 411;
	c->jl3tau  = 0.5*c->jltau;
	c->jlss	= c->hlss;
	c->jl3ss	= c->hl3ss;
	c->ml	  = c->mlss-(c->mlss-c->ml)*exp(-dt/c->mltau);
	c->ml3     = c->ml3ss-(c->ml3ss-c->ml3)*exp(-dt/c->ml3tau);
	c->hl	  = c->hlss-(c->hlss-c->hl)*exp(-dt/c->hltau);
	c->hl3     = c->hl3ss-(c->hl3ss-c->hl3)*exp(-dt/c->hl3tau);
	c->jl	  = c->jlss-(c->jlss-c->jl)*exp(-dt/c->jltau);
	c->jl3     = c->jl3ss-(c->jl3ss-c->jl3)*exp(-dt/c->jl3tau);
	c->inal2   = gnal2*c->ml*c->hl*c->jl*(c->v-c->ena);
	c->inal3   = gnal3*c->ml3*c->hl3*c->jl3*(c->v-c->ena);
	c->inal    = c->inal2 + c->inal3; 
}

void comp_ina (Cell *c)
{
    c->ma	= 0.64*(c->v+37.13)/(1-exp(-0.1*(c->v+37.13)));
	c->mb	= 0.16*exp(-c->v/11);
	if (c->v<-40)
	{
		c->ha = 0.135*exp((70+c->v)/-6.8);
		c->hb = 3.56*exp(0.079*c->v)+310000*exp(0.35*c->v);
		c->ja = (-127140*exp(0.2444*c->v)-0.003474*exp(-0.04391*c->v))*(c->v+37.78)/(1+exp(0.311*(c->v+79.23)));
		c->jb = 0.1212*exp(-0.01052*c->v)/(1+exp(-0.1378*(c->v+40.14)));
	}
	else
	{
		c->ha = 0.0;
		c->hb = 1/(0.13*(1+exp((c->v+10.66)/-11.1)));
		c->ja = 0.0;
		c->jb = 0.3*exp(-0.0000002535*c->v)/(1+exp(-0.1*(c->v+32)));
	}
	c->mtau	= 1/(c->ma+c->mb);
	c->htau	= 1/(c->ha+c->hb);
	c->jtau	= 1/(c->ja+c->jb);
	c->mss	= c->ma*c->mtau;
	c->hss	= c->ha*c->htau;
	c->jss	= 1*c->ja*c->jtau;
	c->m	= c->mss-(c->mss-c->m)*exp(-dt/c->mtau);
	c->h	= c->hss-(c->hss-c->h)*exp(-dt/c->htau);
	c->j	= c->jss-(c->jss-c->j)*exp(-dt/c->jtau);
	c->ina	= gna*pow(c->m,3)*c->h*c->j*(c->v-c->ena);
}

void comp_revs (Cell *c)
{
    c->eca	= (R*temp/(zca*frdy))*log(cao/c->cassl);
	c->ena	= (R*temp/frdy)*log(nao/c->nassl);
	c->ek	= (R*temp/frdy)*log(ko/c->ki);
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
*/