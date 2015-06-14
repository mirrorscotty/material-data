/**
 * @file mechanical.c
 * All sorts of fun stuff related to the Maxwell model for viscoelasticity.
 */

#include "mechanical.h"
#include "pasta.h"
#include "choi-okos.h"
#include "glass-transition.h"
#include <stdlib.h>
#include <math.h>

/**
 * Calculate the shear modulus for a Maxwell material.
 * @param m Set of parameters to use
 * @param t Time [s]
 * @param T Temperature [K]
 * @param M Moisture content [kg/kg db]
 * @returns Shear modulus [Pa]
 */
double MaxwellRelax(maxwell *m, double t, double T, double M)
{
    double E = 0, /* Set the modulus to zero initially */
           tr = ReducedTime(m, t, T, M); /* Reduced time */
    int i; /* Loop index */

    /* Add up each term in the series */
    for(i=0; i<m->n; i++) /* Same problem as in the MeanRelaxTime function */
        E += m->E[i] * exp( -tr/m->tau[i] );

    return E;
}

/**
 * The time derivative of the shear modulus for a Maxwell material.
 * @param m Set of parameters
 * @param t Time [s]
 * @param T Temperature [K]
 * @param M Moisture content [kg/kg db]
 * @returns Shear modulus [Pa]
 */
double DMaxwellRelax(maxwell *m, double t, double T, double M)
{
    double E = 0, /* Set the modulus to zero initially */
           tr = ReducedTime(m, t, T, M); /* Reduced time */
    int i; /* Loop index */

    /* Add up each term in the series */
    for(i=0; i<m->n; i++) /* Same problem as in the MeanRelaxTime function */
        E += -1*m->E[i]/m->tau[i] * exp( -tr/m->tau[i] );

    return E;
}

/**
 * Shear modulus from Rozzi 2002. This version considers the spring constants
 * to be functions of temperature and moisture content. The relaxation times are
 * constant.
 * @param t Time [s]
 * @param T Temperature [K]
 * @param M Moisture content [kg/kg db]
 * @returns Shear modulus [Pa]
 */
double MaxwellRelaxLaura(double t, double T, double M)
{
    double Ea, E1, E2, l1, l2, cg;
    M = M/(1+M); /* Convert from dry basis to wet basis */
    M *= 100; /* The moisture content should be converted to a percentage */

    Ea = 68.18*(1/(1+exp((M-250.92*exp(-0.0091*T))/2.19))+0.078);
    E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238));
    E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    l1 = 7;
    l2 = 110;
    Ea *= 1e6;
    E1 *= 1e6;
    E2 *= 1e6;

    return Ea+E1*exp(-t/l1)+E2*exp(-t/l2);
}

double DMaxwellRelaxLaura(double t, double T, double M)
{
    double Ea, E1, E2, l1, l2, cg;
    M = M/(1+M); /* Convert from dry basis to wet basis */
    M *= 100; /* The moisture content should be converted to a percentage */

    E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238));
    E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    l1 = 7;
    l2 = 110;
    Ea *= 1e6;
    E1 *= 1e6;
    E2 *= 1e6;

    return -E1/l1*exp(-t/l1)-E2/l2*exp(-t/l2);
}

/**
 * Shear modulus from Rozzi 2002. Here, the model is a master curve with both
 * temperature and moisture shift factors.@
 * @param t Time [s]
 * @param T Temperature [K]
 * @param M Moisture content [kg/kg db]
 * @returns Shear modulus [Pa]
 */
double MaxwellRelaxLauraShift(double t, double T, double M)
{
    double aM = exp(104.36*(M+(0.0005*T+0.1569))-10.811),
           aT = exp(0.0806*T-23.849),
           E;
    E = -0.9361+2.51/(1+exp((log(t*aT*aM)-21.01)/5.29));
    return exp(E);
}

