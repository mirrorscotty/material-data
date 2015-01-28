#include "mechanical.h"
#include <math.h>
#include <stdlib.h>

burgers* CreateBurgers()
{
    burgers *b;
    b = (burgers*) calloc(sizeof(burgers), 1);

    b->n = 2;
    b->J = calloc(sizeof(double), b->n);
    b->tau = calloc(sizeof(double), b->n);

    b->J0 = 5.88e-6; // [1/Pa]
    b->J[0] = 1.13e-6; // [1/Pa]
    b->J[1] = 1.02e-6; // [1/Pa]
    b->tau[0] = 2.237; // [s]
    b->tau[1] = 27.485; // [s]
    b->mu0 = 2.96e8; // [Pa s]

    return b;
}

burgers* CreateBurgers1()
{
    burgers *b;
    b = (burgers*) calloc(sizeof(burgers), 1);

    b->n = 2;
    b->J = calloc(sizeof(double), b->n);
    b->tau = calloc(sizeof(double), b->n);

    b->J0 = 5.29e-7; // [1/Pa]
    b->J[0] = 4.86e-8; // [1/Pa]
    b->J[1] = 2.57e-8; // [1/Pa]
    b->tau[0] = 2.397; // [s]
    b->tau[1] = 27.039; // [s]
    b->mu0 = 1.46e10; // [Pa s]

    return b;
}

/**
 * This is a set of data fitted from the parameters in Gina's thesis. Since all
 * of her viscoelasticity parameters are functions of stress, here, they've been
 * fitted to an exponential function.
 * @returns Viscoelasticity parameters
 */
burgerse* CreateBurgersE()
{
    burgerse* b;
    b = (burgerse*) calloc(sizeof(burgerse), 1);

    b->n = 2;
    b->Jb = calloc(sizeof(double), b->n);
    b->Je = calloc(sizeof(double), b->n);
    b->tau = calloc(sizeof(double), b->n);

    b->J0b = .008856; // Fitted to f(x) = b*x^e;
    b->J0e = -.7575;

    b->Jb[0] = .5592; // Fitted to f(x) = b*x^e;
    b->Je[0] = -1.316;

    b->Jb[1] = .1171; // Fitted to f(x) = b*x^e;
    b->Je[1] = -1.189;

    b->tau[0] = 2.286; // Averaged
    b->tau[1] = 27.29;

    b->mu0b = 13150000; // Fitted to f(x) = b*x^e;
    b->mu0e = 1.094;

    return b;
}

void DestroyBurgers(burgers* b)
{
    free(b->J);
    free(b->tau);
    free(b);
}

double BurgersCreep(burgers *b, double t, double T, double M)
{
    double sum = 0;
    int i;

    sum = b->J0 + t/b->mu0;
    for(i=0; i<b->n; i++)
        sum += b->J[i]*(1-exp(-t/b->tau[i]));

    return sum;
}

double DBurgersCreep(burgers *b, double t, double T, double M)
{
    double sum = 0;
    int i;

    sum = 1/b->mu0;
    for(i=0; i<b->n; i++)
        sum += b->J[i]/b->tau[i] * exp(-t/b->tau[i]);

    return sum;
}

/**
 * Creep function for a set of Burgers viscoelasticity parameters. Here, each
 * parameter is a function of applied stress. For stress values less than a
 * somewhat arbitrary value of 10kPa, the pressure is assumed to be equal to
 * 10kPa for purposes of calculating the parameters. This ensures that the
 * function actually returns a reasonable value.
 * @param b Creep function parameters
 * @param t Time [s]
 * @param T Temperature [K]
 * @param M Moisture content [kg/kg db]
 * @param s Applied stress [Pa]
 * @returns Creep function value
 */
double BurgersECreep(burgerse *b, double t, double T, double M, double s)
{
    double sum = 0,
           J0, mu0, *J;
    int i;

    if(s<10000)
        s=10000;

    J = (double*) calloc(sizeof(double), b->n);

    J0 = b->J0b*pow(s, b->J0e);
    mu0 = b->mu0b*pow(s, b->mu0e);
    for(i=0; i<b->n; i++)
        J[i] = b->Jb[i]*pow(s, b->Je[i]);

    sum = J0 + t/mu0;
    for(i=0; i<b->n; i++)
        sum += J[i]*(1-exp(-t/b->tau[i]));

    free(J);

    return sum;
}

/**
 * This is the derivative of the creep function above with respect to time.
 * @param b Creep function parameters
 * @param t Time [s]
 * @param T Temperature [K]
 * @param M Moisture content [kg/kg db]
 * @param s Applied stress [Pa]
 * @returns dJ/dt
 */
double DBurgersECreep(burgerse *b, double t, double T, double M, double s)
{
    double sum = 0,
           mu0, *J;
    int i;

    if(s<10000)
        s=10000;

    J = (double*) calloc(sizeof(double), b->n);

    mu0 = b->mu0b*pow(s, b->mu0e);
    for(i=0; i<b->n; i++)
        J[i] = b->Jb[i]*pow(s, b->Je[i]);

    sum = 1/mu0;
    for(i=0; i<b->n; i++)
        sum += J[i]/b->tau[i] * exp(-t/b->tau[i]);

    free(J);

    return sum;
}

