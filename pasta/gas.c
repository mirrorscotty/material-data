#include <math.h>

#include "pasta.h"

double rho_gas(double wv, double T, double P)
{
    double Mg, /* Molar volume of gas phase */
           R = GASCONST; /* Gas constant */
    return Mg*P/(R*T);
}

/* Vapor pressure of water from 1 C to 374 C (source: Wikipedia)
 * T: Temperature [K]
 * return: vapor pressure [Pa]
 */
double pvap_wat(double T)
{
    double A, B, C, pvap;
    T = T - 273.15;
    
    /* Set Antione equation parameters. The first set is valid from T=1C to 99C,
     * and the second is valid from T=100C to 374C. The two equations are equal
     * at T=108.266 C */
    if(T<108.266) {
        A = 8.07131; B = 1730.63; C = 233.426;
    } else {
        A = 8.14019; B = 1810.94; C = 244.485;
    }

    /* Calcualte vapor pressure with Antione's Equation */
    pvap = pow(10, A - B/(C+T)); 
    pvap = 101300/760; /* Convert from mmHg to Pa */

    return pvap;
}

