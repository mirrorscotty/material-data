#include "material-data.h"
#include "matrix.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    vector *t, *Jlaura0, *Jlaura1, *Jlaura2, *Jlaura3, *Jlaura4;
    matrix *out;
    burgerse *b;
    double ti, Ji;
    int i, n=1000;
    double T=273+40, M=.1, P=.2e6;

    t = linspaceV(0, 1e2, n);
    Jlaura0 = CreateVector(n);
    Jlaura1 = CreateVector(n);
    Jlaura2 = CreateVector(n);
    Jlaura3 = CreateVector(n);
    Jlaura4 = CreateVector(n);

    b = CreateBurgersE();

    for(i=0; i<n; i++) {
        ti = valV(t, i);
        Ji = LLauraCreep(ti, T, M, P);

        Ji = LLauraCreep(ti, T, .1, P);
        setvalV(Jlaura0, i, Ji);
        Ji = LLauraCreep(ti, T, .2, P);
        setvalV(Jlaura1, i, Ji);
        Ji = LLauraCreep(ti, T, .3, P);
        setvalV(Jlaura2, i, Ji);
        Ji = LLauraCreep(ti, T, .4, P);
        setvalV(Jlaura3, i, Ji);
        Ji = LLauraCreep(ti, T, .5, P);
        setvalV(Jlaura4, i, Ji);
    }

    out = CatColVector(6, t, Jlaura0, Jlaura1, Jlaura2, Jlaura3, Jlaura4);

    mtxprntfilehdr(out, "output.csv", "Time\n");

    return;
}


