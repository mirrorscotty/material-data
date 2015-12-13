#include "material-data.h"
#include "matrix.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    vector *e, *xf, *phi;
    matrix *out;
    int i, n=100;
    double ei, xfi, phii, X0 = .33;

    e = linspaceV(0, -.3, n);
    xf = CreateVector(n);
    phi = CreateVector(n);

    for(i=0; i<n; i++) {
        ei = valV(e, i);

        xfi = solidfrac(X0, 300, ei);
        phii = porosity(X0, .15, 300, ei);

        setvalV(xf, i, xfi);
        setvalV(phi, i, phii);
    }

    out = CatColVector(3, e, xf, phi);

    mtxprntfilehdr(out, "output.csv", "strain,solidfrac,porosity\n");

    return;
}


