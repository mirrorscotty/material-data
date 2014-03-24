#include "diff_data.h"
#include <stdlib.h>

diff_data* CreateDiffData()
{
    diff_data *d;
    d = (diff_data*) calloc(sizeof(diff_data), 1);

    d->kw = 1e-18; /* From Zhu et al. 2011 */
    d->muw = .988e-3;
    d->phi = .2;
    d->R = 8.314;
    return d;
}

void DestroyDiffData(diff_data *d)
{
    free(d);
}

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

/* Oswin parameters from Xiong et al. 1991 */
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

/* GAB Data calculated from (Erbas et al. 2005)
 * This data is for semolina.
 */
gab* CreateGABErbas()
{
    gab *d;
    d = (gab*) calloc(sizeof(gab), 1);

    /* Parameters for monolayer moisture content */
    d->m0 = 0.0259;
    d->dHm = 30.9639;
    
    /* C values */
    d->C0 = 170.7282;
    d->dHC = -81.1285;

    /* k values */
    d->k0 = 0.7997;
    d->dHk = -4.1838;

    return d;
}

gab* CreateGABAndrieu()
{
    gab *d;
    d = (gab*) calloc(sizeof(gab), 1);

    /* Parameters for monolayer moisture content */
    d->m0 = 0.0259;
    d->dHm = 30.9639;

    /* C values */
    d->C0 = 170.7282;
    d->dHC = -81.1285;

    /* k values */
    d->k0 = 0.7997;
    d->dHk = -4.1838;

    return d;
}

gab* CreateGABAndrieuK()
{
    gab *d;
    d = (gab*) calloc(sizeof(gab), 1);

    /* Parameters for monolayer moisture content */
    d->m0 = 0.0017;
    d->dHm = 1.1353e3;

    /* C values */
    d->C0 = 2.8339e-8;
    d->dHC = 6.5458e3;

    /* k values */
    d->k0 = 1.8320;
    d->dHk = -238.4933;

    return d;
}


void DestroyGABData(gab *d)
{
    free(d);
}

