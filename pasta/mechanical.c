/**
 * @file mechanical.c
 * All sorts of fun stuff related to the Maxwell model for viscoelasticity.
 */

#include "mechanical.h"
#include "pasta.h"
#include "choi-okos.h"
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
    m->T0 = 298; /* [kg/kg db] */

    /* Moisture shift */
    m->aM0 = -73; /* [s] */
    m->M0 = .14; /* [kg/kg db] */

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
    double aT, aM;

    /* Calculate temperature and moisture shift factors */
    aT = m->aT0*(T - m->T0);
    aM = m->aM0*(M - m->M0);

    return t*aT*aM;
}

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

double MaxwellCreep(maxwell *m, double t, double T, double M)
{
    double J = 0, /* Set the modulus to zero initially */
           tr = ReducedTime(m, t, T, M); /* Reduced time */
    int i; /* Loop index */

    /* Add up each term in the series */
    for(i=0; i<m->n; i++) 
        J += 1/(m->E[i]*m->tau[i])*t + 1/m->E[i];

    return J;
}

/**
 * Derivative of Maxwell material creep compliance function. This isn't
 * actually a function of time, temperature, or moisture; however, the
 * parameters are left there so that it's similar to the normal Maxwell creep
 * function.
 * TODO: Make sure this is mathematically correct with respect to the reduced
 * time thing. It likely isn't.
 */
double DMaxwellCreep(maxwell *m, double t, double T, double M)
{
    double J = 0;
    double i;

    for(i=0; i<m->n; i++)
        J += 1/(m->E[i]*m->tau[i]);

    return J;
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
        stress += MaxwellRelax(m, t-tau, T(tau), M(tau)) * dstrain * dtau;
    }

    return stress;
}

/**
 * Calculate pore pressure based on the Kelvin equation.
 * @param Xdb Moisture content [kg/kg db]
 * @param T Temperature [K]
 * @returns Capillary pressure [Pa]
 */
double pore_press(double Xdb, double T)
{
    double rhow,
           R = GASCONST,
           Vm = 1.802e-5, /* m^3/mol */
           aw;
    choi_okos *co;
    oswin *o;

    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    o = CreateOswinData();
    aw = OswinInverse(o, Xdb, T);
    DestroyOswinData(o);

    return R*T/Vm * log(aw);
}

