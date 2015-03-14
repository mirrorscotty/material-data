/**
 * @file gab.c
 * Implementation of the GAB isotherm equation
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "isotherms.h"
#include "conversions.h"

/**
 * Dealocate a set of GAB parameters
 */
void DestroyGABData(gab *d)
{
    free(d);
}

/**
 * Determine dry basis moisture content based on water activity and temperature
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

