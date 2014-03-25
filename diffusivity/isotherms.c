#include <math.h>

#include "diff_data.h"
#include "isotherms.h"

/* Calculate dry basis moisture content using the Oswin model.
 * aw is water activity
 * T is absolute temperature
 */
double OswinIsotherm(oswin *dat, double aw, double T)
{
    double Xdb;
    T = T-273.15;
    Xdb = (dat->k0 + dat->k1*T) * pow(aw/(1-aw), (dat->n0 + dat->n1*T));
    return Xdb;
}

/* Determine water activity from the Oswin equation.
 * X is dry basis moisture content
 * T is absolute temperature
 */
double OswinInverse(oswin *dat, double X, double T)
{
    double aw;
    T = T-273.15;
    aw = pow(X/(dat->k0 + dat->k1*T), 1/(dat->n0 + dat->n1*T))
        / (pow(X/(dat->k0 + dat->k1*T), 1/(dat->n0 + dat->n1*T)) + 1);
    
    /* Make sure the water activities we're calculating are valid */
    if(aw > 1)
        aw = 1;
    if(aw < 0)
        aw = 0;

    return aw;
}

double OswinDawDx(oswin *dat, double X, double T)
{
    double D, A;
    T = T-273.15;
    A = pow(X/(dat->k0 + dat->k1*T), 1/(dat->n0 + dat->n1*T));
    D = A/(X*(dat->n0 + dat->n1*T)*pow(A+1, 2));
    return D;
}

/* Determine dry basis moisture content based on water activity and temperature
 * using the GAB equation.
 * dat is the set of parameters used to determine the temperature-dependent
 *      constants for the GAB equation.
 * aw is the water activity (between 0 and 1)
 * T is the absolute temperature
 */
double GABIsotherm(gab *dat, double aw, double T)
{
    double Xdb, Xm, C, k;
    //T=T-273.15;
    
    /* Calculate the constant values based on temperature. */
    Xm = dat->m0*exp(dat->dHm/T);
    C = dat->C0*exp(dat->dHC/T);
    k = dat->k0*exp(dat->dHk/T);

    /* Return the moisture content */
    Xdb = Xm * (C*k*aw)/((1-k*aw)*(1-k*aw+C*k*aw));
    return Xdb;
}

/* Determine the water activity based on dry basis moisture content and
 * temperature using the GAB equation.
 * X: dry basis moisture content
 * T: temperature [K]
 */
double GABInverse(gab *dat, double Xdb, double T)
{
    double xm, c, k, A, B, C, y, aw;

    xm = dat->m0*exp(dat->dHm/T);
    c = dat->C0*exp(dat->dHC/T);
    k = dat->k0*exp(dat->dHk/T);

    A = k*(1-c)/(xm*c);
    B = (c-2)/(xm*c);
    C = 1/(xm*c*k);
    y = 1/Xdb;

    B = B-y;

    /* Solve aw/Xdb = A*aw^2 + B*aw + C */
    aw = (-B-sqrt(B*B-4*A*C))/(2*A);

    return aw;
}

