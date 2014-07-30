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
//    sum = 0;
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
