#include "material-data.h"
#include "matrix.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    vector *Xdb, *aw;
    matrix *out;
    oswin *o;
    double awi, Xdbi, T=273.15+80;
    int i, n=100;

    o = CreateOswinData();

    aw = linspaceV(0.05, .95, n);
    Xdb = CreateVector(n);

    for(i=0; i<n; i++) {
        awi = valV(aw, i);
        Xdbi = OswinIsotherm(o, awi, T);

        setvalV(Xdb, i, Xdbi);
    }

    out = CatColVector(2, aw, Xdb);

    mtxprntfilehdr(out, "output.csv", "aw,Xdb\n");

    return;
}


