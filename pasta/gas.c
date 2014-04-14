#include <math.h>

#include "pasta.h"

/** Gas density based on the ideal gas law
 * @param wv Vapor mass fraction [-]
 * @param T Temperature [K]
 * @param P Pressure [Pa]
 * @returns Density [kg/m^3]
 */
double rho_gas(double wv, double T, double P)
{
    double Mg = 1/(MWAIR/(1-wv) + MWWAT/wv), /* Average molar mass of gas */
           R = GASCONST; /* Gas constant */
    return Mg*P/(R*T);
}

/**
 * Calculate the molar density of the gas phase based on the ideal gas law.
 * Source: http://en.wikipedia.org/wiki/Molar_volume
 * @param T Temperature [K]
 * @param P Pressure [Pa]
 * @returns Molar density [mol/m^3]
 */
double moldens_gas(double T, double P)
{
    double R = GASCONST; /* Gas constant */
    return P/(R*T);
}

/**
 * Viscosity of air using Sutherland's formula.
 * This formula is valid as long as the gas is ideal, and the temperature is
 * between 0 and 555K. Pressure should be less than 3.45MPa.
 * Source: http://en.wikipedia.org/wiki/Viscosity
 * @param T Temperature [K]
 * @returns Viscosity [Pa s]
 */
double visc_air(double T)
{
    double T0 = 291.15, /* Reference temperature [K] */
           C = 120, /* Sutherland's constant for air */
           mu0 = 1.827e-5; /* Reference viscosity [Pa s] */

    return mu0*(T0+C)/(T+C) * pow(T/T0, 3/2);
}

/**
 * Vapor pressure of water from 1 C to 374 C (source: Wikipedia)
 * @param T Temperature [K]
 * @returns vapor pressure [Pa]
 */
double pvap_wat(double T)
{
    double A, B, C, pvap;
    T = T - 273.15;
    
    /* Set Antione equation parameters. The first set is valid from T=1C to 99C,
     * and the second is valid from T=100C to 374C. The two equations are equal
     * at T=108.266 C */
    if(T<108.266) {
        A = 8.07131; B = 1730.63; C = 233.426;
    } else {
        A = 8.14019; B = 1810.94; C = 244.485;
    }

    /* Calcualte vapor pressure with Antione's Equation */
    pvap = pow(10, A - B/(C+T)); 
    pvap = pvap * 101300./760.; /* Convert from mmHg to Pa */

    return pvap;
}

