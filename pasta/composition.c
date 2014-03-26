/**
 * @file Equations for determining the composition pasta based on the dependent
 * variables. The independent variables are mass concentration of water, mass
 * fraction of water vapor, porosity, temperature, and pressure.
 */

#include "pasta.h"
#include "choi-okos.h"

/**
 * Calculate total gas phase concentration using gas saturation
 * @param cw mass concentration of water [kg/m^3]
 * @param wv mass fraction of water vapor [-]
 * @param phi porosity
 * @param T temperature [K]
 * @param P pressure [Pa]
 * @returns Gas concentration [kg/m^3]
 */
double conc_gas(double cw, double wv, double phi, double T, double P)
{
    double rhog = rho_gas(wv, T, P); /* Gas density */
    return sat_gas(cw, phi, T) * rhog/phi;
}

/**
 * Vapor mass concentration in the gas phase
 * @param cw mass concentration of water [kg/m^3]
 * @param wv mass fraction of water vapor [-]
 * @param phi porosity [-]
 * @param T temperature [K]
 * @param P pressure [Pa]
 * @returns Vapor concentration [kg/m^3]
 */
double conc_vap(double cw, double wv, double phi, double T, double P)
{
    double cg = conc_gas(cw, wv, phi, T, P); /* Gas concentration */
    return cg*wv;
}

/**
 * Mass concentration of air in the gas phase
 * @param cw mass concentration of water [kg/m^3]
 * @param wv: mass fraction of water vapor [-]
 * @param phi: porosity [-]
 * @param T temperature [K]
 * @param P pressure [Pa]
 * @returns Air concentration [kg/m^3]
 */
double conc_air(double cw, double wv, double phi, double T, double P)
{
    double cg = conc_gas(cw, wv, phi, T, P), /* Gas concentration */
           wa = 1-wv; /* air mass fraction */
    return cg*wa;
}

/**
 * Determine water saturation
 * @param cw mass concentration of water [kg/m^3]
 * @param phi porosity [-]
 * @param T Temperature [K]
 * @returns Saturation [-]
 */
double sat_wat(double cw, double phi, double T)
{
    double rhow; /* Water density */

    /* Calculate the density of water from the Choi-Okos equations */
    choi_okos *co;

    /* Water density */
    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    //return cw*phi/rhow;
    return cw/(phi*rhow);
}

/**
 * Dry basis moisture content
 * @param cw mass concentration water [kg/m^3]
 * @param phi: porosity [-]
 * @returns Moisture content [kg/kg db]
 */
double mdb_wat(double cw, double phi, double T)
{
    double rhos; /* Solid denstiy */
    choi_okos *co;
    
    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    return cw/((1-phi)*rhos);
}

/**
 * Dry basis moisture content when all pores are saturated
 * @param phi porosity [-]
 * @param T temperature [K]
 * @returns Moisture content [kg/kg db]
 */
double mdb_wat_sat(double phi, double T)
{
    double rhos, rhow;
    double cw_sat;
    choi_okos *co;

    /* Water density */
    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    /* Solid density */
    co = CreateChoiOkos(PASTACOMP); 
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    cw_sat = phi*rhow;
    return cw_sat/((1-phi)*rhos);
}

/**
 * Gas saturation
 * @param cw Mass concentration of water [kg/m^3]
 * @param phi Porosity [-]
 * @param T Temperature [K]
 * @return Saturation [-]
 */
double sat_gas(double cw, double phi, double T)
{
    /* The only things filling the pores are water and gas */
    return 1-sat_wat(cw, phi, T);
}

/**
 * Calculate the mole fraction of vapor in the gas phase
 * @param wv Mass fraction of water vapor in gas phase [-]
 * @returns Mole fraction [-]
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

