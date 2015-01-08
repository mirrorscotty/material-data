/**
 * @file poisson-test.c
 * Plots the poisson ratio as a function of porosity.
 */

#include "material-data.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    vector *phi, *v1, *v2, *v3;
    matrix *out;
    double v01=.3, v02=.4, v03=.45,
           phii, v1i, v2i, v3i;
    int i, n=100;


    phi = linspaceV(0, .5, n);
    v1 = CreateVector(n);
    v2 = CreateVector(n);
    v3 = CreateVector(n);

    for(i=0; i<n; i++) {
        phii = valV(phi, i);
        v1i = poisson(phii, v01);
        v2i = poisson(phii, v02);
        v3i = poisson(phii, v03);

        setvalV(v1, i, v1i);
        setvalV(v2, i, v2i);
        setvalV(v3, i, v3i);
    }

    out = CatColVector(4, phi, v1, v2, v3);

    mtxprntfilehdr(out, "output.csv", "phi,v=0.3,v=0.4,v=0.45\n");

    return;
}


