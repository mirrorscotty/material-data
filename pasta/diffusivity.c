#include "../choi-okos/choi-okos.h"
#include "isotherms.h"
#include "constants.h"
#include "conversions.h"
#include "pasta.h"
#include "diffusivity.h"
#include <math.h>

/**
 * Calculate capillary diffusivity in a porous medium. Most of the equations are
 * taken from Zhu 2011.
 * @param d Diffusivity data structure
 * @param o Set of Oswin parameters
 * @param X Moisture content [kg/kg, db]
 * @param T Temperature [K]
 * @returns Diffusivity [m^2/s]
 */
double CapillaryDiff(double X, double T)
{
    double Dcap,
           DPcDSw, 
           DawDX,
           aw, /* Water activity */
           Xs, /* Saturated moisture content */
           kw,
           kwi = KWINTR, /* Intrinsic permeability for water */
           Sw,
           Sr = SR, /* Irreducible water saturation */
           R = GASCONST, /* Gas constant */
           muw = visc_wat(T), /* Water viscosity */
           phi = POROSITY, /* Constant porosity */
           fphi;
    oswin *o;
    choi_okos *co;
    
    o = OSWINDATA();

    /* Water activity from the Oswin isotherm model */
    aw = OswinInverse(o, X, T);

    /* Derivative of water activity with respect to d.b. moisture content */
    DawDX = OswinDawDx(o, X, T);

    /* Estimate Xs by pretending that the material is saturated at aw=.95 */
    Xs = OswinIsotherm(o, .95, T); /* TODO: Fix this */
    //Xs = mdb_wat_sat(phi, T);

    /* Porosity factor */
    fphi = 1;

    /* Water permeability */
    Sw = X/Xs;
    if(Sw>Sr)
        kw = kwi * pow((Sw-Sr)/(1-Sr), 3)*fphi;
    else
        kw = 0;
    
    /* Derivative of capillary pressure with respect to water saturation.
     * Equation from Miranda and Silva 2005 */
    co = CreateChoiOkos(0, 0, 0, 0, 0, 1, 0);
    DPcDSw = -rho(co,T)*R*T/co->MW_wat * 1/aw * Xs * DawDX;
    DestroyChoiOkos(co);

    /* Capillary diffusivity from p104 of Zhu thesis */
    Dcap = -kw/(muw*phi) * DPcDSw;

    return Dcap; 
}

/**
 * Calculate capillary diffusivity in a porous medium. The key equation is taken
 * from Zhu et al. 2011.
 * @param d Diffusivity data structure
 * @param o Set of Oswin parameters
 * @param X Moisture content [kg/kg, db]
 * @param T Temperature [K]
 * @returns Diffusivity [m^2/s]
 */
double CapDiff(double X, double T)
{
    double D, /* Diffusivity */
           DawDe, /* Derivative of water activity with respect to volume frac */
           aw, /* Water activity */
           vl, /* Molar volume of water */
           e, /* Volume fraction of water */
           Xs, /* Maximum dry-basis moisture content */
           Sw, /* Saturation of water */
           Sr = SR, /* Irreducible water saturation */
           kw, /* Permeability of water (kwi*kwr) */
           kwi = KWINTR, /* Intrinsic permeability of water */
           fphi, /* Porosity factor */
           phi = POROSITY, /* Assume constant porosity */
           muw = visc_wat(T), /* Viscosity of water */
           R = GASCONST; /* Gas constant */
    oswin *o;
    choi_okos *co;

    o = OSWINDATA();

    /* Calculate the molar volume of water */
    co = CreateChoiOkos(0, 0, 0, 0, 0, 1, 0);
    vl = co->MW_wat/rho(co, T);
    DestroyChoiOkos(co);

    /* Estimate Xs by pretending that the material is saturated at aw=.95 */
    Xs = OswinIsotherm(o, .95, T);
    //Xs = mdb_wat_sat(phi, T);

    /* Porosity factor */
    fphi = 1;

    /* Water permeability */
    Sw = X/Xs;
    if(Sw>Sr)
        kw = kwi*pow((Sw-Sr)/(1-Sr), 3)*fphi;
    else
        kw = 0;
 
    /* Volume fraction of water */
    e = volfrac_wat(conc_wat(X, POROSITY, T), T);
    e = phi * X/Xs;

    /* Derivative of water activity with respect to volume fraction water */
    DawDe = OswinDawDx(o, X, T) * Xs/phi;
    
    /* Water activity from the Oswin isotherm model */
    aw = OswinInverse(o, X, T);

    /* Equation 3.45 from Zhu et al. 2011 */
    D = e*kw/muw * R*T/vl * (log(aw) + e/aw * DawDe);

    DestroyOswinData(o);

    return D;
}

/**
 * Determine the capillary pressure based on dry basis moisture content and
 * temperature. The equation is from Miranda and Silva 2005
 * @param o set of Oswin isotherm parameters
 * @param X Moisture content [kg/kg, db]
 * @param T Temperature [K]
 * @returns Pressure [Pa]
 */
double CapillaryPressure(oswin *o, double X, double T)
{
    double aw, Pc, R = GASCONST;
    choi_okos *co;
    co = CreateChoiOkos(WATERCOMP);
    /* Calculate water activity from the Oswin isotherm */
    aw = OswinInverse(o, X, T);

    /* Calculate capillary pressure based on water activity */
    Pc = -rho(co,T)*R*T/co->MW_wat * log(aw);
    DestroyChoiOkos(co);
    return Pc;
}

/** 
 * Calculate diffusivity in pasta based on the model outlined in chapter 10 of
 * the Handbook of Food Engineering, Second Edition. This function uses the
 * Oswin isotherm to determine binding energy.
 * @param X Moisture Content [kg/kg db]
 * @param T Temperature [K]
 * @returns Diffusivity [m^2/s]
 *
 * @see DiffCh10GAB
 */
double DiffCh10(double X, double T)
{
    oswin *dat;
    dat = OSWINDATA();

    double Deff,
           D0 = 6.3910e-8, /* [m^2/s] Source: Handbook of Food Engineering */
           //Ea = 25900, /* Source: Litchfield and Okos (1986) */
           Ea = 21760, /* [J/mol] Source: Xiong et al. (1991) */
           K = 1032.558, /* Source: Xiong et al. (1991) */
           Eb = BindingEnergyOswin(dat, X, T),
           R = 8.314; /* Gas Constant */

    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    Deff = D0 * exp(-Ea/(R*T))
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

/**
 * Modification of the diffusivity equation from Xiong et al. 1991 to use
 * the self-diffusivity of water.
 */
double DiffCh10Mod(double X, double T)
{
    oswin *dat;
    dat = OSWINDATA();

    double Deff,
           Dself = SelfDiffWater(T),
           phi = POROSITY,
           tau = 10.5409, /* Fitted value for tortuosity based on 55 C data */
           Sw = sat_wat(conc_wat(X, POROSITY, T), POROSITY, T),
           K = 1032.558, /* Source: Xiong et al. (1991) */
           Eb = BindingEnergyOswin(dat, X, T),
           R = 8.314; /* Gas Constant */

    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    Deff = phi/tau * Dself
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

/** 
 * Calculate diffusivity in pasta based on the model outlined in chapter 10 of
 * the Handbook of Food Engineering, Second Edition. This function uses the
 * GAB isotherm to determine binding energy.
 * @param X Moisture Content [kg/kg db]
 * @param T Temperature [K]
 * @returns Diffusivity [m^2/s]
 *
 * @see DiffCh10
 */
double DiffCh10GAB(double X, double T)
{
    gab *dat;
    dat = GABDATA();

    double Deff,
           D0 = 6.3910e-8, /* Source: Xiong et al (1991) */
           //Ea = 25900, /* Source: Litchfield and Okos (1986) */
           Ea = 21760, /* [J/mol] Source: Xiong et al. (1991) */
           K = 1032.6, /* Source: Xiong et al. (1991) */
           Eb = BindingEnergyGAB(dat, X, T),
           R = 8.314; /* Gas Constant */


    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    Deff = D0 * exp(-Ea/(R*T))
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

/**
 * Calculate the diffusivity of water in a nonpolar gas.
 * Source: Transport Phenomena, Second Edition, Eq 17.2-1
 *
 * Calculates the diffusivity at low pressures, and results in approximately a
 * 6-8% error at atmospheric pressure.
 *
 * @param T Temperature [K]
 * @param P Pressure [Pa]
 * @returns Diffusivity [m^2/s]
 */
double VaporDiff(double T, double P)
{
    double D, /* Binary diffusivity [cm^2/s] */
           Tca = 647.14, /* Critical temperature (Water) [K] */
           Tcb = 126.21, /* Critical temperature (Nitrogen) [K] */
           Pca = 217.8, /* Critical pressure (Water) [atm] */
           Pcb = 33.46, /* Critical pressure (Nitrogen) [atm] */
           Ma = 18.0153, /* Molar mass (Water) [g/mol] */
           Mb = 28.0134, /* Molar mass (Nitrogen) [g/mol] */
           a = 3.64e-4, /* Dimensionless parameter for H2O with a nonpolar gas*/
           b = 2.334; /* Dimensionless parameter for H2O with a nonpolar gas */

    P = PA_TO_ATM(P); /* Convert from Pa to atm */

    /* The temperatures for this equation must be in Kelvin, and the pressures
     * should be in Pascals. The diffusivity is calculated in cm^2/s */
    D = a*pow(T/sqrt(Tca*Tcb), b);
    D *= (pow(Pca*Pcb, 1/3.) * pow(Tca*Tcb, 5/12.) * sqrt(1/Ma+1/Mb))/P;

    D *= 1e-4; /* Convert to m^2/s */

    return D;
}

double VaporDiffCh10(double X, double T)
{
    double D0gas = 2e-5,
           D0liq = 6.3910e-8;
    return D0gas/D0liq * DiffCh10(X, T);
}

/**
 * Temperature-dependent self diffusion coefficent for water.
 * Valid from 0C to 100C. From Holz et al. 2000
 * @param T Temperature [K]
 * @returns Diffusion coefficient [m^2/s]
 */
double SelfDiffWater(double T)
{
    double D,
           D0 = 1.635e-8,
           Ts = 215.05,
           gamma = 2.063;
    D = D0*pow(T/Ts-1, gamma);
    return D;
}

double KnudsenDiff(double T)
{
    double D,
           a = 0, /* Pore radius */
           R = GASCONST, /* Gas constant */
           pi = M_PI, /* Value of pi */
           Mw = 18.0153; /* Molar mass of water (g/mol) */

    D = 8*a/3 * sqrt(R*T/(2*pi*Mw));

    return D;
}
