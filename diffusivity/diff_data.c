#include "diff_data.h"
#include <stdlib.h>

diff_data* CreateDiffData()
{
    diff_data *d;
    d = (diff_data*) calloc(sizeof(diff_data), 1);

    d->kw = 1e-18; /* Intrinsic Permeability From Zhu et al. 2011 */
    d->muw = .988e-3; /* Water viscosity (Zhu et al. 2011) */
    d->phi = .9; /* Regular pasta porosity (Xiong 1991) */
    d->Sr = .02; /* Irreducible water saturation Note: This value is made up */
    d->R = 8.314; /* Gas constant */
    return d;
}

void DestroyDiffData(diff_data *d)
{
    free(d);
}

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

