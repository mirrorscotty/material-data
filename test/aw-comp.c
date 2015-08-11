/**
 * aw-comp.c
 * Compares the Oswin isotherm parameters from Xiong (1992) to the refitted
 * parameters based on her original data.
 */

#include "material-data.h"
#include "matrix.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    vector *Xdb11, *Xdb12, *Xdb21, *Xdb22, *aw;
    matrix *out;
    oswin *o1, *o2;
    double awi, Xdbi, T1=273.15+40, T2=273.15+80;
    int i, n=100;

    o1 = CreateOswinXiong();
    o2 = CreateOswinXiongR();

    aw = linspaceV(0.05, .95, n);
    Xdb11 = CreateVector(n);
    Xdb12 = CreateVector(n);
    Xdb21 = CreateVector(n);
    Xdb22 = CreateVector(n);

    for(i=0; i<n; i++) {
        awi = valV(aw, i);
        Xdbi = OswinIsotherm(o1, awi, T1);
        setvalV(Xdb11, i, Xdbi);
        Xdbi = OswinIsotherm(o1, awi, T2);
        setvalV(Xdb12, i, Xdbi);
        Xdbi = OswinIsotherm(o2, awi, T1);
        setvalV(Xdb21, i, Xdbi);
        Xdbi = OswinIsotherm(o2, awi, T2);
        setvalV(Xdb22, i, Xdbi);
    }

    out = CatColVector(5, aw, Xdb11, Xdb12, Xdb21, Xdb22);

    mtxprntfilehdr(out, "output.csv", "aw,XdbO(40),XdbO(80),XdbR(40),XdbR(80)\n");

    return;
}


