#include "choi-okos.h"
#include "isotherms.h"
#include "constants.h"
#include "conversions.h"
#include "pasta.h"
#include "diffusivity.h"
#include <math.h>

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

