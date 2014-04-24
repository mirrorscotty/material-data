/**
 * @file mechanical.c
 * All sorts of fun stuff related to the Maxwell model for viscoelasticity.
 */

#include "mechanical.h"
#include <stdlib.h>
#include <math.h>

maxwell* CreateMaxwell()
{
    double nterms = 4;
    maxwell *m;

    /* Allocate memory for all the stuff */
    m->E = (double *) calloc(sizeof(double), nterms);
    m->tau = (double *) calloc(sizeof(double), nterms);
    m->n = nterms;

    /* Viscoelastic modulus [Pa] */
    m->E[0] = 6.6e5;
    m->E[1] = 1.8e6;
    m->E[2] = 2.5e6;
    m->E[3] = 1.2e7;

    /* Relaxation time constant [s] */
    m->tau[0] = 4.61e8;
    m->tau[1] = 5.58e7;
    m->tau[2] = 1.32e6;
    m->tau[3] = 1.56e5;

    /* Temperature shift */
    m->aT0 = -0.013; /* [s/K] */
    m->T0 = 298; /* [kg/kg db] */

    /* Moisture shift */
    m->aM0 = -73; /* [s] */
    m->M0 = .14; /* [kg/kg db] */

    return m;
}

void DestroyMaxwell(maxwell *m)
{
    free(m->E);
    free(m->tau);

    free(m);
}

double ReducedTime(maxwell *m, double t, double T, double M)
{
    double aT, aM;

    /* Calculate temperature and moisture shift factors */
    aT = m->aT0*(T - m->T0);
    aM = m->aM0*(M - m->M0);

    return t*aT*aM;
}

double MaxwellModulus(maxwell *m, double t, double T, double M)
{
    double E = 0, /* Set the modulus to zero initially */
           tr = ReducedTime(m, t, T, M); /* Reduced time */
    int i; /* Loop index */

    /* Add up each term in the series */
    for(i=0; i<m->n; i++)
        E += m->E[i] * exp( -tr/m->tau[i] );

    return E;
}

/**
 * Mean relaxation time defined by Vrentas, Jarzebski and Duda (1975).
 * Calculate mean relaxation time using the following formula:
 * \f[
 * \lambda_m = \frac{\int_0^\infty sG(s) ds}{\int_0^\infty G(s) ds}
 * \f]
 * where G(s) is the shear relaxation modulus (Pa) and the denominator is equal
 * to the viscosity at zero shear (Pa s).
 * @param m Set of Maxwell parameters to use.
 * @returns Mean relaxation time [s]
 */
double MeanRelaxTime(maxwell *m)
{
    double G = 0,
           sG = 0;
    int i;

    /* Use the formulas for each infinite integral to evaluate them */
    for(i=0; i<m->n; i++) {
        G += m->tau[i];
        sG += pow(m->tau[i], 2);
    }

    return sG/G;
}

double MaxwellStress(maxwell *m, double t,
        double (*strain)(double), /* Strain vs time */
        double (*M)(double), /* Moisture content vs time [kg/kg db] */
        double (*T)(double)) /* Temperature vs time [K] */
{
    int nsteps = 20, /* Number of steps for numeric integration */
        i; /* loop index */
    double dtau = t/nsteps, /* Calculate the size of each step */
           h = 1e-10, /* For numerical evaluation of derivatives */
           stress = 0, /* Set stress to zero initially */
           dstrain,
           tau;

    /* Integrate using the trapazoid rule from tau = 0 to tau = t */
    for(i=0; i<nsteps; i++) {
        /* Calculate tau from the loop index and the step size */
        tau = i*dtau;
        /* Derivative of strain with respect to time */
        dstrain = (strain(tau+h) - strain(tau-h))/(2*h);
        /* Calculate the stress for this time period and add it in */
        stress += MaxwellModulus(m, t-tau, T(tau), M(tau)) * dstrain * dtau;
    }

    return stress;
}

