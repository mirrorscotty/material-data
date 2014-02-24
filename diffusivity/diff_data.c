#include "diff_data.h"
#include <stdlib.h>

/* Set the parameters for the Oswin isotherm model to the ones determined in
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

void DestroyOswinData(oswin *d)
{
    free(d);
}

/* Set all the parameters for the GAB Equation.
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
    d->m0 = 0.0251;
    d->dHm = 33.1177;
    
    /* C values */
    d->C0 = 3.4334;
    d->dHC = 80.0229;

    /* k values */
    d->k0 = 0.8543;
    d->dHk = -1.9096;

    return d;
}

void DestroyGABData(gab *d)
{
    free(d);
}

