#include "isotherms.h"
#include <stdlib.h>

/**
 * Set of GAB data for potato from Singh 1996. The original data came from
 * Kiranoudis et al. 1993.
 * Valid for 30 C, 45 C, and 60 C
 */
gab* CreateGABPotatoKir()
{
    gab *d;
    double R = 8.314e-3;
    d = (gab*) calloc(sizeof(gab), 1);

    /* Parameters for monolayer moisture content */
    d->m0 = 5.60e-2;
    d->dHm = 1.39/R;
    
    /* C values */
    d->C0 = 0.017;
    d->dHC = 20/R;

    /* k values */
    d->k0 = 0.94;
    d->dHk = -0.92/R;

    return d;
}

/**
 * Potato isotherm data from Chemkhi et al. 2004. It's only valid for
 * temperatures from 40 to 50 C, and does not include temperature effects.
 */
gab* CreateGABPotatoChemkhi()
{
    gab *d;
    d = (gab*) calloc(sizeof(gab), 1);

    /* Parameters for monolayer moisture content */
    d->m0 = 0.209;
    d->dHm = 0.0;
    
    /* C values */
    d->C0 = 4.416;
    d->dHC = 0.0;

    /* k values */
    d->k0 = 0.976;
    d->dHk = 0.0;

    return d;
}

