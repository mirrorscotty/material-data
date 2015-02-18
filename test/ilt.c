#include "matrix.h"
#include "material-data.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

vector* ilt_euler(double complex (*)(double complex), vector*, int);

double complex f(double complex s)
{
    return 1/cpow(s, 2);
}

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

    return sum;
}

double complex LGinaRelax(double complex s)
{
    burgerse  *b;
    double T = 298,
           M = .2,
           P = 10000;
    b = CreateBurgersE();

    return 1/(s*s*LBurgersECreep(b, s, T, M, P));
}

int main(int argc, char *argv[])
{
    vector *t, *result;
    t = CreateVector(5);
    setvalV(t, 0, 1);
    setvalV(t, 1, 10);
    setvalV(t, 2, 100);
    setvalV(t, 3, 1000);
    setvalV(t, 4, 10000);

    result = ilt_euler(&LGinaRelax, t, 32);
    PrintVector(result);

    return 0;
}

