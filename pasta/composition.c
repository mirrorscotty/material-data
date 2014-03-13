#include "pasta.h"
#include "choi-okos.h"

/* Calculate total gas phase concentration using gas saturation
 * cw: mass concentration of water [kg/m^3]
 * wv: mass fraction of water vapor
 * phi: porosity
 * T: temperature [K]
 * P: pressure [Pa]
 */
double conc_gas(double cw, double wv, double phi, double T, double P)
{
    double rhog = rho_gas(wv, T, P);
    return sat_gas(cw, phi, T) * rhog/phi;
}

/* Vapor mass concentration in the gas phase
 * cw: mass concentration of water [kg/m^3]
 * wv: mass fraction of water vapor
 * phi: porosity
 * T: temperature [K]
 * P: pressure [Pa]
 */
double conc_vap(double cw, double wv, double phi, double T, double P)
{
    double cg = conc_gas(cw, wv, phi, T, P); /* Gas concentration */
    return cg*wv;
}

/* Mass concentration of air in the gas phase
 * cw: mass concentration of water [kg/m^3]
 * wv: mass fraction of water vapor
 * phi: porosity
 * T: temperature [K]
 * P: pressure [Pa]
 */
double conc_air(double cw, double wv, double phi, double T, double P)
{
    double cg = conc_gas(cw, wv, phi, T, P), /* Gas concentration */
           wa = 1-wv; /* air mole fraction */
    return cg*wa;
}

/* Determine water saturation
 * cw: mass concentration of water [kg/m^3]
 * phi: porosity
 * T: Temperature [K]
 */
double sat_wat(double cw, double phi, double T)
{
    double rhow; /* Water density */

    /* Calculate the density of water from the Choi-Okos equations */
    choi_okos *co;
    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    return cw*phi/rhow;
}

/* Dry basis moisture content
 * cw: mass concentration water
 * phi: porosity
 */
double mdb_wat(double cw, double phi, double T)
{
    double rhos; /* Solid denstiy */
    choi_okos *co;
    
    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    return cw/(phi*rhos);
}

/* Gas saturation
 * cw: mass concentration of water [kg/m^3]
 * phi: porosity
 * T: Temperature [K]
 */
double sat_gas(double cw, double phi, double T)
{
    /* The only things filling the pores are water and gas */
    return 1-sat_wat(cw, phi, T);
}

/* Calculate the mole fraction of vapor in the gas phase
 * wv is the mass fraction of water vapor in gas phase
 */
double molefrac_vap(double wv)
{
    double Mv = MWWAT, /* Vapor molar mass */
           wa = 1-wv, /* Mass fraction of air */
           Mavg = 1/(MWAIR/wa + MWWAT/wv); /* Average gas phase molar mass */

    /* rhov = rhog * wv, rhov/rhog = wv
     * X = rhov/rhog * Mavg/Mv */
    return wv * Mavg/Mv;
}

