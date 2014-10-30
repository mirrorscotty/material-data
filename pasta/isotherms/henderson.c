/**
 * Implementation of the Henderson isotherm equation as found in
 * (Litchfield 1992)
 */

#include <stdlib.h>
#include <math.h>

#include "isotherms.h"
#include "conversions.h"

/**
 * Set of parameters for the modified Henderson isotherm from Litchfield 1992.
 */
henderson* CreateHendersonData()
{
    henderson *d;
    d = (henderson*) calloc(sizeof(henderson), 1);

    d->A = 1.3638e-11;
    d->B = 2.5728e2;
    d->C = 4.1686;
    d->D = 2.9060;

    return d;
}

void DestroyHendersonData(henderson *d)
{
    free(d);
}

/**
 * Modified Henderson isotherm model from Litchfield 1992.
 * @param d Set of parameters to use
 * @param X Moisture content [kg/kg db]
 * @param T Temperature [K]
 * @returns Water activity [-]
 */
double HendersonIsotherm(henderson *d, double X, double T)
{
    double aw;
    X = X*100;
    aw = 1 - exp(-1*d->A * pow(T - d->B, d->C) * pow(X, d->D));
    return aw;
}

/**
 * Modified Henderson isotherm model from Litchfield 1992. (Inverted)
 * @param d Set of parameters to use
 * @param aw Water activity [-]
 * @param T Temperature [K]
 * @returns Moisture content [kg/kg db]
 */
double HendersonInverse(henderson *d, double aw, double T)
{
    double Xdb;
    Xdb = pow(log(1-aw)/(-1*d->A * pow(T-d->B, d->C)), 1/d->D);
    return Xdb/100;
}

