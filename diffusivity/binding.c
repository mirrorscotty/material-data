#include <math.h>

#include "diff_data.h"
#include "isotherms.h"

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
    Eb = T*T*R*dlnawdT;

    return Eb;
}

/* Calculate binding energy based on the Oswin isotherm model.
 * dat is a set of constants fot the isotherm equation.
 * X is the dry basis moisture content
 * T is absolute temperature
 */
double BindingEnergyOswin(oswin *dat, double X, double T)
{
    double Eb, h, R, dlnawdT;
    /* Set the value of dx used for numerical differentiation */
    h = .000001;
    /* Gas constant */
    R = 8.314;

    dlnawdT = (log(OswinInverse(dat, X, T+h)) - log(OswinInverse(dat, X, T-h)))/(2*h);

    Eb = T*T*R*dlnawdT;

    return Eb;
}

