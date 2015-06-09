#include "isotherms.h"
#include <stdlib.h>

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
 * Parameters fitted from data in Singh et al. 1996
 * Note: This is the same data as from Xiong 1992, but refitted to the GAB
 * isotherm instead of the Oswin isotherm.
 */
gab* CreateGABSingh()
{
    gab *d;
    double R = 8.314e-3;
    d = (gab*) calloc(sizeof(gab), 1);

    /* Parameters for monolayer moisture content */
    d->m0 = 0.0322; //[kg/kg db]
    d->dHm = 1.48/R; //[kJ/mol]

    /* C values */
    d->C0 = 2.370; //[-]
    d->dHC = 14.89/R; //[kJ/mol]

    /* k values */
    d->k0 = 1.38; //[-]
    d->dHk = -1.14/R; //[kJ/mol]

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

