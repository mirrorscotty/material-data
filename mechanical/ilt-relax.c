#include "material-data.h"
#include "matrix.h"
#include "inv-laplace.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

/* Approximate 0 as .001 for derivatives */
#define ZERO 0.001

double complex LBurgersECreep(burgerse *b,
                              double complex s,
                              double T,
                              double M,
                              double P)
{
    double complex sum;
    double J0, mu0, *J;
    int i;

    if(P<10000)
        P=10000;

    J = (double*) calloc(sizeof(double), b->n);

    J0 = b->J0b*pow(P, b->J0e);
    mu0 = b->mu0b*pow(P, b->mu0e);
    for(i=0; i<b->n; i++)
        J[i] = b->Jb[i]*pow(P, b->Je[i]);

    sum = J0/s + 1/(mu0*s*s);
    for(i=0; i<b->n; i++)
        sum += J[i]*(1/s - 1/(s+1/b->tau[i]));

    free(J);

    return sum;
}

double complex _LGinaRelax(double complex s, void *params)
{
    burgerse *b;
    mechdat *d;
    double result, T, M, P;

    d = (mechdat*) params;
    b = CreateBurgersE();

    T = d->T;
    M = d->M;
    P = d->P;

    result = 1/(s*s*LBurgersECreep(b, s, T, M, P));

    DestroyBurgersE(b);

    return result;
}

double complex _DLGinaRelax(double complex s, void *params)
{
    burgerse *b;
    mechdat *d;
    double result, T, M, P;

    d = (mechdat*) params;
    b = CreateBurgersE();

    T = d->T;
    M = d->M;
    P = d->P;

    result = 1/(s*LBurgersECreep(b, s, T, M, P)) - LGinaRelax(ZERO, T, M, P);

    DestroyBurgersE(b);

    return result;
}

double LGinaRelax(double t, double T, double M, double P)
{
    vector *tv, *Gv;
    mechdat *param;
    double G;
    tv = CreateVector(1);
    setvalV(tv, 0, t);

    param = CreateMechDat(T, M, P);

    Gv = ilt_euler(&_LGinaRelax, tv, 32, (void*)param);
    G = valV(Gv, 0);

    DestroyVector(Gv);
    DestroyVector(tv);
    DestroyMechDat(param);

    return G;
}

double DLGinaRelax(double t, double T, double M, double P)
{
    vector *tv, *Gv;
    mechdat *param;
    double G;
    tv = CreateVector(1);
    setvalV(tv, 0, t);

    param = CreateMechDat(T, M, P);

    Gv = ilt_euler(&_DLGinaRelax, tv, 32, (void*)param);
    G = valV(Gv, 0);

    DestroyVector(Gv);
    DestroyVector(tv);
    DestroyMechDat(param);

    return G;
}

