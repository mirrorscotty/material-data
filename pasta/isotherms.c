#include <stdlib.h>
#include <math.h>

#include "isotherms.h"

/**
 * Set the parameters for the Oswin isotherm model to the ones determined in
 * Gina's thesis.
 */
oswin* CreateOswinData()
{
    oswin *d;
    d = (oswin*) calloc(sizeof(oswin), 1);

    d->k0 = 0.1571;
    d->k1 = -0.0012;
    d->n0 = 0.2076;
    d->n1 = 0.0043;

    return d;
}

/**
 * Oswin parameters from Xiong et al. 1991
 */
oswin* CreateOswinXiong()
{
    oswin *d;
    d = (oswin*) calloc(sizeof(oswin), 1);

    d->k0 = 0.176;
    d->k1 = -1.748e-3;
    d->n0 = 0.182;
    d->n1 = 6.946e-3;

    return d;
}

/**
 * Deallocate a set of Oswin parameters
 */
void DestroyOswinData(oswin *d)
{
    free(d);
}

/**
 * Set all the parameters for the GAB Equation.
 * Data fitted from Gina's thesis.
 * Each set of values describes the temperature dependence of one of the GAB
 * parameters via an equation of the form m = m0*exp(dHm/T), where T is in
 * degrees Celcius.
 */
gab* CreateGABData()
{
    gab *d;
    d = (gab*) calloc(sizeof(gab), 1);

    /* Parameters for monolayer moisture content */
    d->m0 = 3.8004e-4; //0.0251;
    d->dHm = 1.6074e3; //33.1177;
    
    /* C values */
    d->C0 = 4.1929e-4; //3.4334;
    d->dHC = 3.5155e3; //80.0229;

    /* k values */
    d->k0 = 1.0013; //0.8543;
    d->dHk = -65.3646; //-1.9096;

    return d;
}

/**
 * GAB Data calculated from (Erbas et al. 2005)
 * This data is for semolina.
 */
gab* CreateGABErbas()
{
    gab *d;
    d = (gab*) calloc(sizeof(gab), 1);

    /* Parameters for monolayer moisture content */
    d->m0 = 1.3841e-5; //0.0259;
    d->dHm = 2.6428e3; //30.9639;
    
    /* C values */
    d->C0 = 1.2762e12; //170.7282;
    d->dHC = -7.8642e3; //-81.1285;

    /* k values */
    d->k0 = 2.1566; //0.7997;
    d->dHk = -348.9801; //-4.1838;

    return d;
}

/**
 * Parameters fitted from data in Andrieu et al. 1985
 * Actual data copied from Gina's thesis
 */
gab* CreateGABAndrieu()
{
    gab *d;
    d = (gab*) calloc(sizeof(gab), 1);

    /* Parameters for monolayer moisture content */
    d->m0 = 0.0017; //0.0259;
    d->dHm = 1.1353e3; //30.9639;

    /* C values */
    d->C0 = 2.8339e-8; //170.7282;
    d->dHC = 6.5458e3; //-81.1285;

    /* k values */
    d->k0 = 1.8320; //0.7997;
    d->dHk = -238.4933; //-4.1838;

    return d;
}

/**
 * Dealocate a set of GAB parameters
 */
void DestroyGABData(gab *d)
{
    free(d);
}

/**
 * Calculate dry basis moisture content using the Oswin model.
 * @param dat Set of Oswin parameters
 * @param aw Water activity [-]
 * @param T absolute temperature [K]
 * @returns Moisture content [kg/kg db]
 *
 * @see OswinInverse
 */
double OswinIsotherm(oswin *dat, double aw, double T)
{
    double Xdb;

    T = T-273.15; /* Convert from Kelvin to Celcius */

    /* Source: Handbook of Food Engineering, Second Edition (Ch10, Eq 9) */
    Xdb = (dat->k0 + dat->k1*T) * pow(aw/(1-aw), (dat->n0 + dat->n1*T));

    return Xdb;
}

/**
 * Determine water activity from the Oswin equation.
 * @param dat Set of Oswin parameters
 * @param X Dry basis moisture content [kg/kg db]
 * @param T Absolute temperature [K]
 * @returns Water activity [-]
 */
double OswinInverse(oswin *dat, double X, double T)
{
    double aw;
    T = T-273.15; /* Convert from Kelvin to Celcius */

    /* Solved the Oswin equation for water activity using Wolfram|Alpha */
    aw = pow(X/(dat->k0 + dat->k1*T), 1/(dat->n0 + dat->n1*T))
        / (pow(X/(dat->k0 + dat->k1*T), 1/(dat->n0 + dat->n1*T)) + 1);
    
    /* Make sure the water activities we're calculating are valid */
    if(aw > 1)
        aw = 1;
    if(aw < 0)
        aw = 0;

    return aw;
}

/**
 * Calculate the deriviative of water activity with respect to moisture content
 * using the Oswin model. Used for diffusivity calculations.
 * @param dat Set of Oswin parameters
 * @param X Moisture content [kg/kg db]
 * @param T Temperature [K]
 * @returns Deriviative of water activity w.r.t. moisture content
 *
 * @see OswinInverse
 */
double OswinDawDx(oswin *dat, double X, double T)
{
    double D, A;
    T = T-273.15; /* Convert from K to C */

    /* Derivative taken using Wolfram|Alpha */
    A = pow(X/(dat->k0 + dat->k1*T), 1/(dat->n0 + dat->n1*T));
    D = A/(X*(dat->n0 + dat->n1*T)*pow(A+1, 2));

    return D;
}

/** Determine dry basis moisture content based on water activity and temperature
 * using the GAB equation.
 * @param dat Set of parameters used to determine the temperature dependent
 *      constants for the GAB equation.
 * @param aw water activity (between 0 and 1) [-]
 * @param T Absolute temperature [K]
 * @returns Moisture content [kg/kg db]
 *
 * @see GABInverse
 */
double GABIsotherm(gab *dat, double aw, double T)
{
    double Xdb, Xm, C, k;
    
    /* Calculate the constant values based on temperature. Here, temperature
     * should be measured in Kelvins */
    Xm = dat->m0*exp(dat->dHm/T); /* Monolayer moisture content */
    C = dat->C0*exp(dat->dHC/T); /* Guggenheim constant */
    k = dat->k0*exp(dat->dHk/T); /* Factor correcting properties of multilayer
                                    molecules with respect to the bulk liquid */

    /* Return the moisture content */
    Xdb = Xm * (C*k*aw)/((1-k*aw)*(1-k*aw+C*k*aw));

    return Xdb;
}

/**
 * Determine the water activity based on dry basis moisture content and
 * temperature using the GAB equation.
 * @param X dry basis moisture content [kg/kg db]
 * @param T temperature [K]
 * @returns Water activity [-]
 *
 * @see GABIsotherm
 */
double GABInverse(gab *dat, double Xdb, double T)
{
    double xm, c, k, A, B, C, y, aw;

    /* Calculate the constant values based on temperature. */
    xm = dat->m0*exp(dat->dHm/T);
    c = dat->C0*exp(dat->dHC/T);
    k = dat->k0*exp(dat->dHk/T);

    /* Since it's easiest to calculate water activity when the equation is in
     * a quadratic form, determine the constants in front of each power of aw.
     */
    A = k*(1-c)/(xm*c);
    B = (c-2)/(xm*c);
    C = 1/(xm*c*k);

    /* This is really aw/Xdb, but since it's just going to be subtracted from
     * B*aw anyway, the water activity is factored out. */
    y = 1/Xdb;
    B = B-y;

    /* Solve aw/Xdb = A*aw^2 + B*aw + C using the quadratic equation */
    aw = (-B-sqrt(B*B-4*A*C))/(2*A);

    return aw;
}

