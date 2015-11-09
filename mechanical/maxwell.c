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
 * Create a set of maxwell parameters. The data here comes from Cummings et al.
 * 1993 and is for extruded durum semolina pasta.
 * @returns Set of Maxwell material parameters for extruded durum semolina
 */
maxwell* CreateMaxwell()
{
    int nterms = 4;
    maxwell *m;
    m = (maxwell*) calloc(sizeof(maxwell), 1);

    /* Allocate memory for all the stuff */
    m->E = (double *) calloc(sizeof(double), nterms);
    m->tau = (double *) calloc(sizeof(double), nterms);

    /* Set the number of maxwell elements we'll be using */
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
    m->T0 = 298; /* [K] */

    /* Moisture shift */
    m->aM0 = -73; /* [s] */
    m->M0 = .14; /* [kg/kg db] */

    return m;
}

/**
 * @returns Set of Maxwell material parameters for extruded durum semolina
 */
maxwell* CreateMaxwellZhu()
{
    int nterms = 6;
    maxwell *m;
    m = (maxwell*) calloc(sizeof(maxwell), 1);

    /* Allocate memory for all the stuff */
    m->E = (double *) calloc(sizeof(double), nterms);
    m->tau = (double *) calloc(sizeof(double), nterms);

    /* Set the number of maxwell elements we'll be using */
    m->n = nterms;

    /* Viscoelastic modulus [Pa] */
    m->E[0] = 0.6709e6;
    m->E[1] = 0.5312e6;
    m->E[2] = 0.3694e6;
    m->E[3] = 0.3216e6;
    m->E[4] = 0.2613e6;
    m->E[5] = 1.1322e6;

    /* Relaxation time constant [s] */
    m->tau[0] = 0.03;
    m->tau[1] = 0.422;
    m->tau[2] = 3.9;
    m->tau[3] = 98;
    m->tau[4] = 2000;
    m->tau[5] = 299990;

    /* Temperature shift */
    m->aT0 = 0; /* [s/K] */
    m->T0 = 298; /* [K] */

    /* Moisture shift */
    m->aM0 = 9.7577; /* [s] */
    m->M0 = -1.8638; /* [kg/kg db] */

    return m;
}

/**
 * Deallocate a data structure containing maxwell parameters
 * @param m Data structure to delete
 */
void DestroyMaxwell(maxwell *m)
{
    free(m->E);
    free(m->tau);

    free(m);
}

double TimeShift(maxwell *m, double T, double M)
{
    double aT, aM;

    /* Calculate temperature and moisture shift factors */
    aT = exp(m->aT0*(T - m->T0));
    aM = exp(m->aM0*(M - m->M0));

    return aT*aM;
}

/**
 * Calculate the reduced time based on the temperature and moisture shift
 * factors. This adjusts the viscoelastic modulus to account for changes in
 * temperature and moisture content.
 * @param m Set of Maxwell parameters to use
 * @param t Normal (non-reduced) time) [s]
 * @param T Temperature [K]
 * @param M Moisture content [kg/kg db]
 * @returns Reduced time [s]
 */
double ReducedTime(maxwell *m, double t, double T, double M)
{
    return t*TimeShift(m, T, M);
}

/**
 * Mean relaxation time defined by Vrentas, Jarzebski and Duda (1975).
 * Calculate mean relaxation time using the following formula:
 * \f[
 * \lambda_m = \frac{\int_0^\infty sG(s) ds}{\int_0^\infty G(s) ds}
 * \f]
 * where G(s) is the shear relaxation modulus (Pa) and the denominator is equal
 * to the viscosity at zero shear (Pa s).
 *
 * TODO: For some reason, m->n isn't being read properly, and always comes out
 * as zero. Fix this!
 *
 * @param m Set of Maxwell parameters to use.
 * @returns Mean relaxation time [s]
 */
double MeanRelaxTime(maxwell *m)
{
    double G = 0,
           sG = 0;
    int i;

    /* Use the formulas for each infinite integral to evaluate them */
    /* Replace 4 by m->n */
    for(i=0; i<m->n; i++) {
        G += m->tau[i];
        sG += pow(m->tau[i], 2);
    }

    return sG/G;
}

