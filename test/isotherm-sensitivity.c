#include "material-data.h"
#include "matrix.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    vector *aw, *Xdb, *Xdbk0m, *Xdbk0p, *Xdbk1m, *Xdbk1p, *Xdbn0m, *Xdbn0p, *Xdbn1m, *Xdbn1p;
    matrix *out;
    oswin *o;
    double awi, Xdbi, T;
    int i, n=100;

    if(argc != 2) {
        printf("Usage:\n"
               "isotherm-sensitivity <T>\n"
               "<T>\tTemperature (K)\n");
        exit(0);
    }

    T = atoi(argv[1]);

    o = CreateOswinData();

    aw = linspaceV(0.05, .95, n);
    Xdb = CreateVector(n);
    Xdbk0m = CreateVector(n);
    Xdbk0p = CreateVector(n);
    Xdbk1m = CreateVector(n);
    Xdbk1p = CreateVector(n);
    Xdbn0m = CreateVector(n);
    Xdbn0p = CreateVector(n);
    Xdbn1m = CreateVector(n);
    Xdbn1p = CreateVector(n);

    for(i=0; i<n; i++) {
        awi = valV(aw, i);

        o = CreateOswinData();
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdb, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k0 = o->k0*.9;
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbk0m, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k0 = o->k0*1.1;
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbk0p, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k1 = o->k1*.9;
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbk1m, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k1 = o->k1*1.1;
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbk1p, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->n0 = o->n0*.9;
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbn0m, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->n0 = o->n0*1.1;
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbn0p, i, Xdbi);
        DestroyOswinData(o);
        
        o = CreateOswinData();
        o->n1 = o->n1*.9;
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbn1m, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->n1 = o->n1*1.1;
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbn1p, i, Xdbi);
        DestroyOswinData(o);
    }

    out = CatColVector(10, aw, Xdb, Xdbk0m, Xdbk0p, Xdbk1m, Xdbk1p, Xdbn0m, Xdbn0p, Xdbn1m, Xdbn1p);

    mtxprntfilehdr(out, "output.csv", "aw,Xdb,k0-10%,k0+10%,k1-10%,k1+10%,n0-10%,n0+10%,n1-10%,n1+10%\n");

    return;
}

