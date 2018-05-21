#include "../include/luorudy.h"

void LuoRudy::initConst (double CONSTANTS[], double RATES[], double STATES[])
{
    STATES[0] = -84.3801107371;
    STATES[1] = 0.00171338077730188;
    STATES[2] = 0.982660523699656;
    STATES[3] = 0.989108212766685;
    STATES[4] = 0.00017948816388306;
    STATES[5] = 0.00302126301779861;
    STATES[6] = 0.999967936476325;
    STATES[7] = 0.0417603108167287;
    CONSTANTS[0] = 8314;
    CONSTANTS[1] = 310;
    CONSTANTS[2] = 96484.6;
    CONSTANTS[3] = 1;
    CONSTANTS[4] = 100;
    CONSTANTS[5] = 9000;
    CONSTANTS[6] = 1000;
    CONSTANTS[7] = 2;
    CONSTANTS[8] = -25.5;
    CONSTANTS[9] = 23;
    CONSTANTS[10] = 140;
    CONSTANTS[11] = 18;
    CONSTANTS[12] = 0.01833;
    CONSTANTS[13] = 5.4;
    CONSTANTS[14] = 145;
    CONSTANTS[15] = 0.0183;
    CONSTANTS[16] = -59.87;
    CONSTANTS[17] = 0.03921;
    CONSTANTS[18] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[10]/CONSTANTS[11]);
    CONSTANTS[19] =  0.282000* pow((CONSTANTS[13]/5.40000), 1.0 / 2);
    CONSTANTS[20] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log((CONSTANTS[13]+ CONSTANTS[12]*CONSTANTS[10])/(CONSTANTS[14]+ CONSTANTS[12]*CONSTANTS[11]));
    CONSTANTS[21] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[13]/CONSTANTS[14]);
    CONSTANTS[22] =  0.604700* pow((CONSTANTS[13]/5.40000), 1.0 / 2);
    CONSTANTS[23] = CONSTANTS[21];
}

void LuoRudy::compVariables (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{
    ALGEBRAIC[1] = ( 0.320000*(STATES[0]+47.1300))/(1.00000 - exp( - 0.100000*(STATES[0]+47.1300)));
    ALGEBRAIC[8] =  0.0800000*exp(- STATES[0]/11.0000);
    ALGEBRAIC[2] = (STATES[0]<- 40.0000 ?  0.135000*exp((80.0000+STATES[0])/- 6.80000) : 0.00000);
    ALGEBRAIC[9] = (STATES[0]<- 40.0000 ?  3.56000*exp( 0.0790000*STATES[0])+ 310000.*exp( 0.350000*STATES[0]) : 1.00000/( 0.130000*(1.00000+exp((STATES[0]+10.6600)/- 11.1000))));
    ALGEBRAIC[3] = (STATES[0]<- 40.0000 ? ( ( - 127140.*exp( 0.244400*STATES[0]) -  3.47400e-05*exp( - 0.0439100*STATES[0]))*(STATES[0]+37.7800))/(1.00000+exp( 0.311000*(STATES[0]+79.2300))) : 0.00000);
    ALGEBRAIC[10] = (STATES[0]<- 40.0000 ? ( 0.121200*exp( - 0.0105200*STATES[0]))/(1.00000+exp( - 0.137800*(STATES[0]+40.1400))) : ( 0.300000*exp( - 2.53500e-07*STATES[0]))/(1.00000+exp( - 0.100000*(STATES[0]+32.0000))));
    ALGEBRAIC[4] = ( 0.0950000*exp( - 0.0100000*(STATES[0] - 5.00000)))/(1.00000+exp( - 0.0720000*(STATES[0] - 5.00000)));
    ALGEBRAIC[11] = ( 0.0700000*exp( - 0.0170000*(STATES[0]+44.0000)))/(1.00000+exp( 0.0500000*(STATES[0]+44.0000)));
    ALGEBRAIC[5] = ( 0.0120000*exp( - 0.00800000*(STATES[0]+28.0000)))/(1.00000+exp( 0.150000*(STATES[0]+28.0000)));
    ALGEBRAIC[12] = ( 0.00650000*exp( - 0.0200000*(STATES[0]+30.0000)))/(1.00000+exp( - 0.200000*(STATES[0]+30.0000)));
    ALGEBRAIC[6] = ( 0.000500000*exp( 0.0830000*(STATES[0]+50.0000)))/(1.00000+exp( 0.0570000*(STATES[0]+50.0000)));
    ALGEBRAIC[13] = ( 0.00130000*exp( - 0.0600000*(STATES[0]+20.0000)))/(1.00000+exp( - 0.0400000*(STATES[0]+20.0000)));
    ALGEBRAIC[14] = 7.70000 -  13.0287*log(STATES[4]/1.00000);
    ALGEBRAIC[15] =  0.0900000*STATES[5]*STATES[6]*(STATES[0] - ALGEBRAIC[14]);
    ALGEBRAIC[0] = (t>=CONSTANTS[4]&&t<=CONSTANTS[5]&&(t - CONSTANTS[4]) -  floor((t - CONSTANTS[4])/CONSTANTS[6])*CONSTANTS[6]<=CONSTANTS[7] ? CONSTANTS[8] : 0.00000);
    ALGEBRAIC[7] =  CONSTANTS[9]*pow(STATES[1], 3.00000)*STATES[2]*STATES[3]*(STATES[0] - CONSTANTS[18]);
    ALGEBRAIC[16] = (STATES[0]>- 100.000 ? ( 2.83700*(exp( 0.0400000*(STATES[0]+77.0000)) - 1.00000))/( (STATES[0]+77.0000)*exp( 0.0400000*(STATES[0]+35.0000))) : 1.00000);
    ALGEBRAIC[17] =  CONSTANTS[19]*STATES[7]*ALGEBRAIC[16]*(STATES[0] - CONSTANTS[20]);
    ALGEBRAIC[18] = 1.02000/(1.00000+exp( 0.238500*((STATES[0] - CONSTANTS[21]) - 59.2150)));
    ALGEBRAIC[19] = ( 0.491240*exp( 0.0803200*((STATES[0]+5.47600) - CONSTANTS[21]))+ 1.00000*exp( 0.0617500*(STATES[0] - (CONSTANTS[21]+594.310))))/(1.00000+exp( - 0.514300*((STATES[0] - CONSTANTS[21])+4.75300)));
    ALGEBRAIC[20] = ALGEBRAIC[18]/(ALGEBRAIC[18]+ALGEBRAIC[19]);
    ALGEBRAIC[21] =  CONSTANTS[22]*ALGEBRAIC[20]*(STATES[0] - CONSTANTS[21]);
    ALGEBRAIC[22] = 1.00000/(1.00000+exp((7.48800 - STATES[0])/5.98000));
    ALGEBRAIC[23] =  CONSTANTS[15]*ALGEBRAIC[22]*(STATES[0] - CONSTANTS[23]);
    ALGEBRAIC[24] =  CONSTANTS[17]*(STATES[0] - CONSTANTS[16]);
}

void LuoRudy::compRates (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{
    RATES[0] =  (- 1.00000/CONSTANTS[3])*(ALGEBRAIC[0]+ALGEBRAIC[7]+ALGEBRAIC[15]+ALGEBRAIC[17]+ALGEBRAIC[21]+ALGEBRAIC[23]+ALGEBRAIC[24]);
    RATES[1] =  ALGEBRAIC[1]*(1.00000 - STATES[1]) -  ALGEBRAIC[8]*STATES[1];
    RATES[2] =  ALGEBRAIC[2]*(1.00000 - STATES[2]) -  ALGEBRAIC[9]*STATES[2];
    RATES[3] =  ALGEBRAIC[3]*(1.00000 - STATES[3]) -  ALGEBRAIC[10]*STATES[3];
    RATES[4] =  (- 0.000100000/1.00000)*ALGEBRAIC[15]+ 0.0700000*(0.000100000 - STATES[4]);
    RATES[5] =  ALGEBRAIC[4]*(1.00000 - STATES[5]) -  ALGEBRAIC[11]*STATES[5];
    RATES[6] =  ALGEBRAIC[5]*(1.00000 - STATES[6]) -  ALGEBRAIC[12]*STATES[6];
    RATES[7] =  ALGEBRAIC[6]*(1.00000 - STATES[7]) -  ALGEBRAIC[13]*STATES[7];
}

void LuoRudy::compRates_FE (double t, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{

}

void LuoRudy::RL (double DT, double CONSTANTS[], double RATES[], double STATES[], double ALGEBRAIC[])
{

}

/*
   There are a total of 25 entries in the algebraic variable array.
   There are a total of 8 entries in each of the rate and state variable arrays.
   There are a total of 24 entries in the constant variable array.
 */

/*
 * t is time in component environment (millisecond).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is R in component membrane (joule_per_kilomole_kelvin).
 * CONSTANTS[1] is T in component membrane (kelvin).
 * CONSTANTS[2] is F in component membrane (coulomb_per_mole).
 * CONSTANTS[3] is C in component membrane (microF_per_cm2).
 * ALGEBRAIC[0] is I_stim in component membrane (microA_per_cm2).
 * ALGEBRAIC[7] is i_Na in component fast_sodium_current (microA_per_cm2).
 * ALGEBRAIC[15] is i_si in component slow_inward_current (microA_per_cm2).
 * ALGEBRAIC[17] is i_K in component time_dependent_potassium_current (microA_per_cm2).
 * ALGEBRAIC[21] is i_K1 in component time_independent_potassium_current (microA_per_cm2).
 * ALGEBRAIC[23] is i_Kp in component plateau_potassium_current (microA_per_cm2).
 * ALGEBRAIC[24] is i_b in component background_current (microA_per_cm2).
 * CONSTANTS[4] is stim_start in component membrane (millisecond).
 * CONSTANTS[5] is stim_end in component membrane (millisecond).
 * CONSTANTS[6] is stim_period in component membrane (millisecond).
 * CONSTANTS[7] is stim_duration in component membrane (millisecond).
 * CONSTANTS[8] is stim_amplitude in component membrane (microA_per_cm2).
 * CONSTANTS[9] is g_Na in component fast_sodium_current (milliS_per_cm2).
 * CONSTANTS[18] is E_Na in component fast_sodium_current (millivolt).
 * CONSTANTS[10] is Nao in component ionic_concentrations (millimolar).
 * CONSTANTS[11] is Nai in component ionic_concentrations (millimolar).
 * STATES[1] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[2] is h in component fast_sodium_current_h_gate (dimensionless).
 * STATES[3] is j in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[1] is alpha_m in component fast_sodium_current_m_gate (per_millisecond).
 * ALGEBRAIC[8] is beta_m in component fast_sodium_current_m_gate (per_millisecond).
 * ALGEBRAIC[2] is alpha_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[9] is beta_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[3] is alpha_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[10] is beta_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[14] is E_si in component slow_inward_current (millivolt).
 * STATES[4] is Cai in component intracellular_calcium_concentration (millimolar).
 * STATES[5] is d in component slow_inward_current_d_gate (dimensionless).
 * STATES[6] is f in component slow_inward_current_f_gate (dimensionless).
 * ALGEBRAIC[4] is alpha_d in component slow_inward_current_d_gate (per_millisecond).
 * ALGEBRAIC[11] is beta_d in component slow_inward_current_d_gate (per_millisecond).
 * ALGEBRAIC[5] is alpha_f in component slow_inward_current_f_gate (per_millisecond).
 * ALGEBRAIC[12] is beta_f in component slow_inward_current_f_gate (per_millisecond).
 * CONSTANTS[19] is g_K in component time_dependent_potassium_current (milliS_per_cm2).
 * CONSTANTS[20] is E_K in component time_dependent_potassium_current (millivolt).
 * CONSTANTS[12] is PR_NaK in component time_dependent_potassium_current (dimensionless).
 * CONSTANTS[13] is Ko in component ionic_concentrations (millimolar).
 * CONSTANTS[14] is Ki in component ionic_concentrations (millimolar).
 * STATES[7] is X in component time_dependent_potassium_current_X_gate (dimensionless).
 * ALGEBRAIC[16] is Xi in component time_dependent_potassium_current_Xi_gate (dimensionless).
 * ALGEBRAIC[6] is alpha_X in component time_dependent_potassium_current_X_gate (per_millisecond).
 * ALGEBRAIC[13] is beta_X in component time_dependent_potassium_current_X_gate (per_millisecond).
 * CONSTANTS[21] is E_K1 in component time_independent_potassium_current (millivolt).
 * CONSTANTS[22] is g_K1 in component time_independent_potassium_current (milliS_per_cm2).
 * ALGEBRAIC[20] is K1_infinity in component time_independent_potassium_current_K1_gate (dimensionless).
 * ALGEBRAIC[18] is alpha_K1 in component time_independent_potassium_current_K1_gate (per_millisecond).
 * ALGEBRAIC[19] is beta_K1 in component time_independent_potassium_current_K1_gate (per_millisecond).
 * CONSTANTS[23] is E_Kp in component plateau_potassium_current (millivolt).
 * CONSTANTS[15] is g_Kp in component plateau_potassium_current (milliS_per_cm2).
 * ALGEBRAIC[22] is Kp in component plateau_potassium_current (dimensionless).
 * CONSTANTS[16] is E_b in component background_current (millivolt).
 * CONSTANTS[17] is g_b in component background_current (milliS_per_cm2).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[1] is d/dt m in component fast_sodium_current_m_gate (dimensionless).
 * RATES[2] is d/dt h in component fast_sodium_current_h_gate (dimensionless).
 * RATES[3] is d/dt j in component fast_sodium_current_j_gate (dimensionless).
 * RATES[5] is d/dt d in component slow_inward_current_d_gate (dimensionless).
 * RATES[6] is d/dt f in component slow_inward_current_f_gate (dimensionless).
 * RATES[7] is d/dt X in component time_dependent_potassium_current_X_gate (dimensionless).
 * RATES[4] is d/dt Cai in component intracellular_calcium_concentration (millimolar).
 */