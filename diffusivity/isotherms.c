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

/* Calculate binding energy based on the Oswin isotherm model.
 * dat is a set of constants fot the isotherm equation.
 * X is the dry basis moisture content
 * T is absolute temperature
 */
double BindingEnergyOswin(oswin *dat, double X, double T)
{
    double Eb, h, R, dlnawdT;
    h = .00001;
    R = 8.314;

    dlnawdT = (log(OswinInverse(dat, X, T+h)) - log(OswinInverse(dat, X, T-h)))/(2*h);

    Eb = T*T*R*dlnawdT;

    return Eb;
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
/* Wrong */
double GABInverse3(gab *dat, double X, double T)
{
    double xm, c, k, A, B, C, y, aw;
    T = T-273.15;

    xm = dat->m0*exp(dat->dHm/T);
    c = dat->m0*exp(dat->dHC/T);
    k = dat->k0*exp(dat->dHk/T);

    A = 1/(k*c);
    B = (k-2)/k;
    C = (c-k*c)/k;
    y = X/xm;

    aw = (-1*sqrt(-4*A*C*y*y + B*B*y*y - 2*B*y + 1) - B*y + 1)/(2*C*y);

    return aw;
}

/* Also wrong */
double GABInverse2(gab *dat, double X, double T)
{
    double Xm, C, K, aw;
    T = T-273.15;

    Xm = dat->m0*exp(dat->dHm/T);
    C = dat->m0*exp(dat->dHC/T);
    K = dat->k0*exp(dat->dHk/T);

    aw = -1*(C*(Xm-X)+sqrt(C)*sqrt(C*pow(Xm-X, 2)+4*Xm*X) + 2*X)/(2*(C-1)*K*X);
    return aw;
}

/* Determine the water activity based on dry basis moisture content and
 * temperature using the GAB equation.
 * X is dry basis moisture content
 * T is absolute temperature
 */
double GABInverse(gab *dat, double X, double T)
{
    double Xm, C, K, aw;
    T = T-273.15;

    Xm = dat->m0*exp(dat->dHm/T);
    C = dat->m0*exp(dat->dHC/T);
    K = dat->k0*exp(dat->dHk/T);

    aw = (C*(X-Xm)+sqrt(C)*sqrt(C*pow(Xm-X, 2)+4*Xm*X) - 2*X)/(2*(C-1)*K*X);
    return aw;
}

/* Calculate the binding energy based on the GAB equation.
 * X is the dry-basis moisture content
 * T is the absolute temperature
 */
double BindingEnergyGAB(gab *dat, double X, double T)
{
    double Eb, h, R, dlnawdT;
    /* Set the value of dx used for numerical differentiation */
    h = .000001;
    /* Gas constant */
    R = 8.314;

    dlnawdT = (log(GABInverse(dat, X, T+h)) - log(GABInverse(dat, X, T-h)))/(2*h);
    Eb = R*T*T*dlnawdT;

    return Eb;
}

double BindingEnergyGABInt(gab *dat, double X, double T)
{
    double Eb, R, T2, aw1, aw2, h;
    /* Set the value of dx used for numerical differentiation */
    h = 1;
    /* Gas constant */
    R = 8.314;
    T2 = T + h;

    aw1 = GABInverse(dat, X, T);
    aw2 = GABInverse(dat, X, T2);
    
    Eb = R*log(aw2/aw1)/(1/T-1/T2);

    return Eb;
}

