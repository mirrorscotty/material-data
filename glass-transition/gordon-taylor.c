/**
 * @file gordon-taylor.c
 * Functions to calculate glass transition temperature using the Gordon-Taylor
 * equation.
 */

#include "glass-transition.h"
#include <stdlib.h>

gordontaylor* GTSemolina()
{
    gordontaylor *gt;
    gt = (gordontaylor*) calloc(sizeof(gordontaylor), 1);

    gt->Tg1 = 435; /* [K] Tg of solid at zero moisture content */
    gt->Tg2 = 138; /* [K] Tg of water */
    gt->kGT = 3.4; /* Gordon-Taylor constant */

    return gt;
}

void DestroyGT(gordontaylor *gt)
{
    free(gt);
    return;
}

/**
 * Gordon-Taylor equation for predicting glass transition temperature.
 * @param gt Equation parameters
 * @param Xdb Moisture content [kg/kg db]
 * @returns Glass transition temperature [K]
 */
double GordonTaylor(gordontaylor *gt, double Xdb)
{
    double w1, w2, Tg;
    w2 = Xdb/(1+Xdb); /* Convert from dry basis to wet basis moisture content */
    w1 = 1-w2; /* Calculate the mass fraction of solids */

    /* Gordon-Taylor equation */
    Tg = (w1*gt->Tg1 + gt->kGT*w2*gt->Tg2)/(w1+gt->kGT*w2);
    return Tg;
}

/**
 * Use the gordon taylor equation to predict moisture content from glass
 * transition temperature.
 * @param gt Equation parameters
 * @param T Glass transition temperature [K]
 * @returns Moisture content [kg/kg db]
 */
double GordonTaylorInv(gordontaylor *gt, double T)
{
    double Xwb;
    /* Gordon-Taylor equation inverted to calculate wet-basis moisture
     * content. */
    Xwb = (gt->Tg1-T)/(gt->Tg1-gt->Tg2*gt->kGT+(gt->kGT-1)*T);
    /* Convert to dry basis */
    return Xwb/(1-Xwb);
}

