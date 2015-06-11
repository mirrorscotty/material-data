#include "material-data.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[])
{
    vector *Ga, *G1, *G2, *Xdb;
    matrix *out;
    double Xdbi, T=313, Gi;
    int i, n=100;

    Xdb = linspaceV(0.05, .5, n);
    Ga = CreateVector(n);
    G1 = CreateVector(n);
    G2 = CreateVector(n);

    for(i=0; i<n; i++) {
        Xdbi = valV(Xdb, i);
        Xdbi = Xdbi/(Xdbi+1); /* Convert to wet basis */
        Xdbi = Xdbi * 100; /* Convert to a percentage */
        Gi = 68.18*(1/(1+exp((Xdbi-250.92*exp(-0.0091*T))/2.19))+0.078) * 1e6;
        setvalV(Ga, i, Gi);
        Gi = 20.26*exp(-0.0802*(Xdbi+0.0474*T-14.238)) * 1e6;
        setvalV(G1, i, Gi);
        Gi = 2.484 + 6.576/(1+exp((Xdbi-19.36)/0.848)) * 1e6;
        setvalV(G2, i, Gi);
    }

    out = CatColVector(4, Xdb, Ga, G1, G2);

    mtxprntfilehdr(out, "plot-rozzi-xdb.csv", "Xdb,Ga,G1,G2\n");

    return;
}

