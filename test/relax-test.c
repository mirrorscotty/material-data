#include "material-data.h"
#include "matrix.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

vector* ilt_euler(double complex (*)(double complex), vector*, int);

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
    double T = 273+40,
           M = .08,
           P = 200000; //pore_press(M, T);
    b = CreateBurgersE();

    return 1/(s*s*LBurgersECreep(b, s, T, M, P));
}

int main(int argc, char *argv[])
{
    vector *t, *Gcummings, *Glaura, *Ggina;
    matrix *out;
    maxwell *m;
    double ti, Gcummingsi, Glaurai;
    int i, n=1000;
    double T=298, M=.1;

    m = CreateMaxwell();

    t = linspaceV(0, 1e8, n);
    Gcummings = CreateVector(n);
    Glaura = CreateVector(n);

    for(i=0; i<n; i++) {
        ti = valV(t, i);
        Gcummingsi = MaxwellRelax(m, ti, T, M);
        Glaurai = MaxwellRelaxLaura(ti, T, M);

        setvalV(Gcummings, i, Gcummingsi);
        setvalV(Glaura, i, Glaurai);
    }

    Ggina = ilt_euler(&LGinaRelax, t, 32);

    out = CatColVector(4, t, Gcummings, Glaura, Ggina);

    mtxprntfilehdr(out, "output.csv", "Time,Cummings,Rozzi,Bressani\n"); 
    return;
}


