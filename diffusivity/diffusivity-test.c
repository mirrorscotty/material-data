#include "diffusivity.h"
/**
 * @file diffusivity.c
 * Several models for calculating diffusivity in porous media.
 */

#include "choi-okos.h"
#include "isotherms.h"
#include "constants.h"
#include "conversions.h"
#include "pasta.h"
#include "diffusivity.h"
#include "binding.h"
#include <math.h>
#include <stdlib.h>

double DiffCh10_test(double T)
{

    double Deff, R = GASCONST;
    DiffXiongData *d;
    d = CreateDefaultXiongData();

    Deff = d->D0 * exp(-d->Ea/(R*T));
    DestroyXiongData(d);

    return Deff;
}

/**
 * Modification of the diffusivity equation from Xiong et al. 1991 to use
 * the self-diffusivity of water.
 */
double DiffCh10mod_test(double T)
{
    oswin *dat;
    dat = OSWINDATA();

    double Deff,
           Dself = SelfDiffWater(T),
           phi = POROSITY,
           tau = 10.5409, /* Fitted value for tortuosity based on 55 C data */
           R = GASCONST;

    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    Deff = phi/tau * Dself;

    return Deff;
}
