#include "material-data.h"
#include "matrix.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    vector *e, *xf;
    matrix *out;
    int i, n=100;
    double ei, xfi;

    e = linspaceV(0, -.3, n);
    xf = CreateVector(n);

    for(i=0; i<n; i++) {
        ei = valV(e, i);
        xfi = solidfrac(.33, 300, ei);

        setvalV(xf, i, xfi);
    }

    out = CatColVector(2, e, sf);

    mtxprntfilehdr(out, "output.csv", "strain,solidfrac\n");

    return;
}


