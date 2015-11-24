#include "material-data.h"
#include "matrix.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    vector *e, *sf;
    matrix *out;
    int i, n=100;
    double ei, sfi;

    e = linspaceV(0, -.3, n);
    sf = CreateVector(n);

    for(i=0; i<n; i++) {
        ei = valV(e, i);
        sfi = solidfrac(.33, 300, ei);

        setvalV(sf, i, sfi);
    }

    out = CatColVector(2, e, sf);

    mtxprntfilehdr(out, "output.csv", "strain,solidfrac\n");

    return;
}


