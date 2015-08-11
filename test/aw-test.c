/**
 * aw-test.c
 * Print a table of isotherm data to "output.csv"
 */

#include "material-data.h"
#include "matrix.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    vector *Xdb1, *Xdb2, *aw;
    matrix *out;
    oswin *o;
    double awi, Xdbi, T1=273.15+40, T2=273.15+80;
    int i, n=100;

    o = CreateOswinData();

    aw = linspaceV(0.05, .95, n);
    Xdb1 = CreateVector(n);
    Xdb2 = CreateVector(n);

    for(i=0; i<n; i++) {
        awi = valV(aw, i);
        Xdbi = OswinIsotherm(o, awi, T1);
        setvalV(Xdb1, i, Xdbi);
        Xdbi = OswinIsotherm(o, awi, T2);
        setvalV(Xdb2, i, Xdbi);
    }

    out = CatColVector(3, aw, Xdb1, Xdb2);

    mtxprntfilehdr(out, "output.csv", "aw,Xdb(40),Xdb(80)\n");

    return;
}


