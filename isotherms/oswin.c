/**
 * @file oswin.c
 * Implementation of the Oswin isotherm equation
 */
#include <stdlib.h>
#include <math.h>

#include "isotherms.h"
#include "conversions.h"

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

