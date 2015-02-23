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

double MaxwellCreepLaura(double t, double T, double M)
{
    double Ea, E1, E2, l1, l2, A;

    Ea = 68.18*(1/(1+exp((M-250.92*exp(-0.0091*T))/2.19))+0.078);
    E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238));
    E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    l1 = 7;
    l2 = 110;
    Ea *= 1e6;
    E1 *= 1e6;
    E2 *= 1e6;

    A = (E1+Ea)*l1+(E2+Ea)*l2;
    return 1/Ea * (1-exp(-Ea*t/A));
}

double DMaxwellCreepLaura(double t, double T, double M)
{
    double Ea, E1, E2, l1, l2, A;

    Ea = 68.18*(1/(1+exp((M-250.92*exp(-0.0091*T))/2.19))+0.078);
    E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238));
    E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    l1 = 7;
    l2 = 110;
    Ea *= 1e6;
    E1 *= 1e6;
    E2 *= 1e6;

    A = (E1+Ea)*l1+(E2+Ea)*l2;
    return 1/A * exp(-Ea*t/A);
}

/**
 * Use stress relaxation data to create a creep compliance function. This was
 * solved for using Maple for a two-element Maxwell solid. The values for the
 * viscoelastic parameters are calculated using the equations from Laura's
 * thesis.
 */
double MaxwellCreepConverted(double t, double T, double M)
{
    double Ea, E1, E2, l1, l2, cg;

    Ea = 68.18*(1/(1+exp((M-250.92*exp(-0.0091*T))/2.19))+0.078);
    E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238));
    E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    l1 = 7;
    l2 = 110;
    Ea *= 1e6;
    E1 *= 1e6;
    E2 *= 1e6;

    /* Solved using Maple */
    cg = (-l2 * E1 * E2 + l2 * E1 * Ea - E1 * l1 * (E1 + E2 + Ea) - l2 * E2 *
            E2 + l1 * E2 * Ea - l2 * E2 * Ea) * exp(-(E1 * l1 + E2 * l2 + l1 *
                    Ea + l2 * Ea) / l1 / l2 / (E1 + E2 + Ea) * t / 0.2e1) *
            sinh(t / l1 / l2 / (E1 + E2 + Ea) * sqrt(E1 * E1 * l1 * l1 + 0.2e1
                        * E1 * E2 * l1 * l2 + 0.2e1 * E1 * Ea * l1 * l1 - 0.2e1
                        * E1 * Ea * l1 * l2 + l2 * l2 * E2 * E2 - 0.2e1 * E2 *
                        Ea * l1 * l2 + 0.2e1 * E2 * Ea * l2 * l2 + Ea * Ea * l1
                        * l1 - 0.2e1 * Ea * Ea * l1 * l2 + Ea * Ea * l2 * l2) /
                    0.2e1) / (E1 + E2 + Ea) * pow(E1 * E1 * l1 * l1 + 0.2e1 *
                        E1 * E2 * l1 * l2 + 0.2e1 * E1 * Ea * l1 * l1 - 0.2e1 *
                        E1 * Ea * l1 * l2 + l2 * l2 * E2 * E2 - 0.2e1 * E2 * Ea
                        * l1 * l2 + 0.2e1 * E2 * Ea * l2 * l2 + Ea * Ea * l1 *
                        l1 - 0.2e1 * Ea * Ea * l1 * l2 + Ea * Ea * l2 * l2,
                        -0.1e1 / 0.2e1) / Ea + (-(E1 + E2) * cosh(t / l1 / l2 /
                                (E1 + E2 + Ea) * sqrt(E1 * E1 * l1 * l1 + 0.2e1
                                    * E1 * E2 * l1 * l2 + 0.2e1 * E1 * Ea * l1
                                    * l1 - 0.2e1 * E1 * Ea * l1 * l2 + l2 * l2
                                    * E2 * E2 - 0.2e1 * E2 * Ea * l1 * l2 +
                                    0.2e1 * E2 * Ea * l2 * l2 + Ea * Ea * l1 *
                                    l1 - 0.2e1 * Ea * Ea * l1 * l2 + Ea * Ea *
                                    l2 * l2) / 0.2e1) * exp(-(E1 * l1 + E2 * l2
                                        + l1 * Ea + l2 * Ea) / l1 / l2 / (E1 +
                                            E2 + Ea) * t / 0.2e1) + E1 + E2 +
                            Ea) / (E1 + E2 + Ea) / Ea;

    return cg;
}

double DMaxwellCreepConverted(double t, double T, double M)
{
    double Ea, E1, E2, l1, l2, cg;

    Ea = 68.18*(1/(1+exp((M-250.92*exp(-0.0091*T))/2.19))+0.078);
    E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238));
    E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    l1 = 7;
    l2 = 110;
    Ea *= 1e6;
    E1 *= 1e6;
    E2 *= 1e6;


    cg = -(-E1 * l2 * E2 + E1 * l2 * Ea - E1 * l1 * (E1 + E2 + Ea) - l2 * E2 * E2 + l1 * E2 * Ea - l2 * E2 * Ea) * (E1 * l1 + E2 * l2 + l1 * Ea + l2 * Ea) / l1 / l2 * pow(E1 + E2 + Ea, -0.2e1) * exp(-(E1 * l1 + E2 * l2 + l1 * Ea + l2 * Ea) / l1 / l2 / (E1 + E2 + Ea) * t / 0.2e1) * sinh(t / l1 / l2 / (E1 + E2 + Ea) * sqrt(E1 * E1 * l1 * l1 + 0.2e1 * E1 * E2 * l1 * l2 + 0.2e1 * E1 * Ea * l1 * l1 - 0.2e1 * E1 * Ea * l1 * l2 + l2 * l2 * E2 * E2 - 0.2e1 * E2 * Ea * l1 * l2 + 0.2e1 * E2 * Ea * l2 * l2 + Ea * Ea * l1 * l1 - 0.2e1 * Ea * Ea * l1 * l2 + Ea * Ea * l2 * l2) / 0.2e1) * pow(E1 * E1 * l1 * l1 + 0.2e1 * E1 * E2 * l1 * l2 + 0.2e1 * E1 * Ea * l1 * l1 - 0.2e1 * E1 * Ea * l1 * l2 + l2 * l2 * E2 * E2 - 0.2e1 * E2 * Ea * l1 * l2 + 0.2e1 * E2 * Ea * l2 * l2 + Ea * Ea * l1 * l1 - 0.2e1 * Ea * Ea * l1 * l2 + Ea * Ea * l2 * l2, -0.1e1 / 0.2e1) / Ea / 0.2e1 + (-E1 * l2 * E2 + E1 * l2 * Ea - E1 * l1 * (E1 + E2 + Ea) - l2 * E2 * E2 + l1 * E2 * Ea - l2 * E2 * Ea) * exp(-(E1 * l1 + E2 * l2 + l1 * Ea + l2 * Ea) / l1 / l2 / (E1 + E2 + Ea) * t / 0.2e1) * cosh(t / l1 / l2 / (E1 + E2 + Ea) * sqrt(E1 * E1 * l1 * l1 + 0.2e1 * E1 * E2 * l1 * l2 + 0.2e1 * E1 * Ea * l1 * l1 - 0.2e1 * E1 * Ea * l1 * l2 + l2 * l2 * E2 * E2 - 0.2e1 * E2 * Ea * l1 * l2 + 0.2e1 * E2 * Ea * l2 * l2 + Ea * Ea * l1 * l1 - 0.2e1 * Ea * Ea * l1 * l2 + Ea * Ea * l2 * l2) / 0.2e1) / l1 / l2 * pow(E1 + E2 + Ea, -0.2e1) / Ea / 0.2e1 + (-(E1 + E2) * sinh(t / l1 / l2 / (E1 + E2 + Ea) * sqrt(E1 * E1 * l1 * l1 + 0.2e1 * E1 * E2 * l1 * l2 + 0.2e1 * E1 * Ea * l1 * l1 - 0.2e1 * E1 * Ea * l1 * l2 + l2 * l2 * E2 * E2 - 0.2e1 * E2 * Ea * l1 * l2 + 0.2e1 * E2 * Ea * l2 * l2 + Ea * Ea * l1 * l1 - 0.2e1 * Ea * Ea * l1 * l2 + Ea * Ea * l2 * l2) / 0.2e1) / l1 / l2 / (E1 + E2 + Ea) * sqrt(E1 * E1 * l1 * l1 + 0.2e1 * E1 * E2 * l1 * l2 + 0.2e1 * E1 * Ea * l1 * l1 - 0.2e1 * E1 * Ea * l1 * l2 + l2 * l2 * E2 * E2 - 0.2e1 * E2 * Ea * l1 * l2 + 0.2e1 * E2 * Ea * l2 * l2 + Ea * Ea * l1 * l1 - 0.2e1 * Ea * Ea * l1 * l2 + Ea * Ea * l2 * l2) * exp(-(E1 * l1 + E2 * l2 + l1 * Ea + l2 * Ea) / l1 / l2 / (E1 + E2 + Ea) * t / 0.2e1) / 0.2e1 + (E1 + E2) * cosh(t / l1 / l2 / (E1 + E2 + Ea) * sqrt(E1 * E1 * l1 * l1 + 0.2e1 * E1 * E2 * l1 * l2 + 0.2e1 * E1 * Ea * l1 * l1 - 0.2e1 * E1 * Ea * l1 * l2 + l2 * l2 * E2 * E2 - 0.2e1 * E2 * Ea * l1 * l2 + 0.2e1 * E2 * Ea * l2 * l2 + Ea * Ea * l1 * l1 - 0.2e1 * Ea * Ea * l1 * l2 + Ea * Ea * l2 * l2) / 0.2e1) * (E1 * l1 + E2 * l2 + l1 * Ea + l2 * Ea) / l1 / l2 / (E1 + E2 + Ea) * exp(-(E1 * l1 + E2 * l2 + l1 * Ea + l2 * Ea) / l1 / l2 / (E1 + E2 + Ea) * t / 0.2e1) / 0.2e1) / (E1 + E2 + Ea) / Ea;

    return cg;
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
    int i;

    for(i=0; i<m->n; i++)
        J += 1/(m->E[i]*m->tau[i]);

    return J;
}

