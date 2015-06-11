#include "material-data.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[])
{
    vector *Ga, *G1, *G2, *T;
    matrix *out;
    double X=.1, Ti, Gi;
    int i, n=100;

    T = linspaceV(298, 353, n);
    Ga = CreateVector(n);
    G1 = CreateVector(n);
    G2 = CreateVector(n);

    X = X/(1+X); /* Convert to wet basis */
    X = X * 100; /* Convert to a percentage */
    for(i=0; i<n; i++) {
        Ti = valV(T, i);

        Gi = 68.18*(1/(1+exp((X-250.92*exp(-0.0091*Ti))/2.19))+0.078) * 1e6;
        setvalV(Ga, i, Gi);
        Gi = 20.26*exp(-0.0802*(X+0.0474*Ti-14.238)) * 1e6;
        setvalV(G1, i, Gi);
        Gi = 2.484 + 6.576/(1+exp((X-19.36)/0.848)) * 1e6;
        setvalV(G2, i, Gi);
    }

    out = CatColVector(4, T, Ga, G1, G2);

    mtxprntfilehdr(out, "plot-rozzi-t.csv", "T,Ga,G1,G2\n");

    return;
}

