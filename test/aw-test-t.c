/**
 * @file aw-test-t.c
 * Plot water activity vs. temperature.
 */

#include "material-data.h"
#include "matrix.h"
#include <stdio.h>

void PlotOswin()
{
    vector *aw1, *aw2, *aw3, *aw4, *aw5, *T;
    matrix *out;
    oswin *o;
    double Xdb1=.05, Xdb2=.10, Xdb3=.15, Xdb4=.20, Xdb5 = .25, Ti, awi;
    int i, n=100;

    o = CreateOswinData();

    T = linspaceV(298, 353, n);
    aw1 = CreateVector(n);
    aw2 = CreateVector(n);
    aw3 = CreateVector(n);
    aw4 = CreateVector(n);
    aw5 = CreateVector(n);

    for(i=0; i<n; i++) {
        Ti = valV(T, i);
        awi = OswinInverse(o, Xdb1, Ti);
        setvalV(aw1, i, awi);
        awi = OswinInverse(o, Xdb2, Ti);
        setvalV(aw2, i, awi);
        awi = OswinInverse(o, Xdb3, Ti);
        setvalV(aw3, i, awi);
        awi = OswinInverse(o, Xdb4, Ti);
        setvalV(aw4, i, awi);
        awi = OswinInverse(o, Xdb5, Ti);
        setvalV(aw5, i, awi);
    }

    out = CatColVector(6, T, aw1, aw2, aw3, aw4, aw5);

    mtxprntfilehdr(out, "oswin-t.csv", "T(K),X=.05,X=.10,X=.15,X=.20,X=.25\n");

    return;
}

void PlotGAB()
{
    vector *aw1, *aw2, *aw3, *aw4, *aw5, *T;
    matrix *out;
    gab *o;
    double Xdb1=.05, Xdb2=.10, Xdb3=.15, Xdb4=.20, Xdb5 = .25, Ti, awi;
    int i, n=100;

    o = CreateGABData();

    T = linspaceV(298, 353, n);
    aw1 = CreateVector(n);
    aw2 = CreateVector(n);
    aw3 = CreateVector(n);
    aw4 = CreateVector(n);
    aw5 = CreateVector(n);

    for(i=0; i<n; i++) {
        Ti = valV(T, i);
        awi = GABInverse(o, Xdb1, Ti);
        setvalV(aw1, i, awi);
        awi = GABInverse(o, Xdb2, Ti);
        setvalV(aw2, i, awi);
        awi = GABInverse(o, Xdb3, Ti);
        setvalV(aw3, i, awi);
        awi = GABInverse(o, Xdb4, Ti);
        setvalV(aw4, i, awi);
        awi = GABInverse(o, Xdb5, Ti);
        setvalV(aw5, i, awi);
    }

    out = CatColVector(6, T, aw1, aw2, aw3, aw4, aw5);

    mtxprntfilehdr(out, "gab-t.csv", "T(K),X=.05,X=.10,X=.15,X=.20,X=.25\n");

    return;
}

int main(int argc, char *argv[])
{
    PlotGAB();
    PlotOswin();
    return;
}

