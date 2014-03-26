#include "pasta.h"
#include "isotherms.h"
#include "diff_data.h"

/**
 * Calculate the rate of evaporation
 * @param cw Mass concentration of water [kg/m^3]
 * @param phi Porosity [-]
 * @param T temperature [K]
 * @returns Evaporation rate
 */
double evap(double cw, double wv, double phi, double T, double P)
{
    double I,
           kevap = 1, /* Evaporation rate constant [1/s] */
           R = GASCONST, /* Gas constant */
           Mv = P/R*T, /* Molar volume of vapor (from ideal gas law) */
           pveq, /* Equilibrium vapor pressure */
           pv = P * molefrac_vap(wv), /* Vapor partial pressure */
           Xdb; /* Dry basis moisture content */
    oswin *o;

    Xdb = mdb_wat(cw, phi, T);

    o = CreateOswinData();
    pveq = OswinInverse(o, Xdb, T) * pvap_wat(T);
    DestroyOswinData(o);

    I = kevap * Mv/(R*T) * (pveq - pv);
    return I;
}

