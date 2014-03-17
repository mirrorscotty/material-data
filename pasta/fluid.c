#include <math.h>

#include "pasta.h"
#include "isotherms.h"
#include "diff_data.h"

/* Viscosity of water (from 0 C to 370 C)
 * T is temperature in units of Kelvins
 * mu is viscosity with units of N*s/m^2
 */
double visc_wat(double T)
{
    double mu;
    mu = 2.414e-5 * pow(10, 247.8/(T-140));
    return mu;
}

double perm_wat(double cw, double phi, double T)
{
    double kw, /* Permeability of water */
           kwi = PERMWAT, /* Intrinsic permeability of water */
           fphi, /* Porosity factor */
           Xs, /* Saturated moisture content (d.b.) */
           Sw, /* Water saturation */
           Sr; /* Irreducible water saturation */
    oswin *o;
    o = CreateOswinData();

    fphi = 1;
    Xs = mdb_wat_sat(phi, T); /* Xdb when pores are saturated with water */
    Sr = 0.05/Xs; /* Completely made up number! */
    Sw = sat_wat(cw, phi, T);

    DestroyOswinData(o);

    if(Sw>Sr)
        kw = kwi*pow((Sw-Sr)/(1-Sr), 3) * fphi;
    else
        kw = 0;
    return kw;
}

