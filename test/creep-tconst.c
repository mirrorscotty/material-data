/**
 * @creep-tconst.c
 * Plot the stress relaxation data from Rozzi (2002) as a function of moisture
 * content at constant times. This is then compared against the creep function
 * (calculated numerically) based on the same data.
 */

#include "material-data.h"
#include "matrix.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

double Err(double Ji, double Gi)
{
    return fabs(Gi-1/Ji)/Gi;
}

int main(int argc, char *argv[])
{
    vector *M, *Jlaura0, *Jlaura1, *Jlaura2, *Jlaura3, *Jlaura4;
    vector *Glaura0, *Glaura1, *Glaura2, *Glaura3, *Glaura4;
    vector *Err0, *Err1, *Err2, *Err3, *Err4;
    matrix *out;
    burgerse *b;
    double ti, Ji, Gi;
    int i, n=1000;
    double Ti=273+40, t=20, P=.2e6, Mi;

    M = linspaceV(0.005, .5, n);
    Jlaura0 = CreateVector(n);
    Jlaura1 = CreateVector(n);
    Jlaura2 = CreateVector(n);
    Jlaura3 = CreateVector(n);
    Jlaura4 = CreateVector(n);

    Glaura0 = CreateVector(n);
    Glaura1 = CreateVector(n);
    Glaura2 = CreateVector(n);
    Glaura3 = CreateVector(n);
    Glaura4 = CreateVector(n);

    Err0 = CreateVector(n);
    Err1 = CreateVector(n);
    Err2 = CreateVector(n);
    Err3 = CreateVector(n);
    Err4 = CreateVector(n);

    for(i=0; i<n; i++) {
        Mi = valV(M, i);

        Ji = LLauraCreep(1*t, Ti, Mi, P);
        setvalV(Jlaura0, i, Ji);
        Ji = LLauraCreep(2*t, Ti, Mi, P);
        setvalV(Jlaura1, i, Ji);
        Ji = LLauraCreep(3*t, Ti, Mi, P);
        setvalV(Jlaura2, i, Ji);
        Ji = LLauraCreep(4*t, Ti, Mi, P);
        setvalV(Jlaura3, i, Ji);
        Ji = LLauraCreep(5*t, Ti, Mi, P);
        setvalV(Jlaura4, i, Ji);

        Gi = MaxwellRelaxLaura(1*t, Ti, Mi);
        setvalV(Glaura0, i, Gi);
        Gi = MaxwellRelaxLaura(2*t, Ti, Mi);
        setvalV(Glaura1, i, Gi);
        Gi = MaxwellRelaxLaura(3*t, Ti, Mi);
        setvalV(Glaura2, i, Gi);
        Gi = MaxwellRelaxLaura(4*t, Ti, Mi);
        setvalV(Glaura3, i, Gi);
        Gi = MaxwellRelaxLaura(5*t, Ti, Mi);
        setvalV(Glaura4, i, Gi);

        setvalV(Err0, i, Err(valV(Jlaura0, i), valV(Glaura0, i)));
        setvalV(Err1, i, Err(valV(Jlaura1, i), valV(Glaura1, i)));
        setvalV(Err2, i, Err(valV(Jlaura2, i), valV(Glaura2, i)));
        setvalV(Err3, i, Err(valV(Jlaura3, i), valV(Glaura3, i)));
        setvalV(Err4, i, Err(valV(Jlaura4, i), valV(Glaura4, i)));

    }

    out = CatColVector(16, M, Jlaura0, Jlaura1, Jlaura2, Jlaura3, Jlaura4, Glaura0, Glaura1, Glaura2, Glaura3, Glaura4, Err0, Err1, Err2, Err3, Err4);

    mtxprntfilehdr(out, "output.csv", "Time,J(t=20s),J(t=40s),J(t=60s),J(t=80s),J(t=100s),G(t=20s),G(t=40s),G(t=60s),G(t=80s),G(t=100s),Err(t=20s),Err(t=40s),Err(t=60s),Err(t=80s),Err(t=100s)\n");

    return;
}

