#include "material-data.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double l2errornorm(vector *x, vector *f1, vector *f2)
{
    int i;
    double norm = 0, base = 0,
           dx, fa, fb;
    for(i=0; i<len(x-1); i++) {
        dx = valV(x, i+1) - valV(x, i);
        fa = pow(valV(f1, i) - valV(f2, i), 2);
        fb = pow(valV(f1, i+1) - valV(f2, i+1), 2);
        norm += dx*(fb+fa);

        fa = pow(valV(f1, i), 2);
        fb = pow(valV(f1, i+1), 2);
        base += dx*(fb+fa);
    }
    return sqrt((.5*norm)/(.5*base));
}

int main(int argc, char *argv[])
{
    vector *aw, *Xdb, *XdbTm, *XdbTp, *Xdbk0m, *Xdbk0p, *Xdbk1m, *Xdbk1p, *Xdbn0m, *Xdbn0p, *Xdbn1m, *Xdbn1p;
    matrix *out;
    oswin *o;
    double awi, Xdbi, T, dev;
    int i, n=100;
    char *hdr;

    if(argc != 3) {
        printf("Usage:\n"
               "isotherm-sensitivity <T> <dev>\n"
               "<T>\tTemperature (K)\n"
               "<dev>\tAmount to change each input variable by\n");
        exit(0);
    }

    T = atof(argv[1]);
    dev = atof(argv[2]);

    o = CreateOswinData();

    aw = linspaceV(0.005, .95, n);
    Xdb = CreateVector(n);
    XdbTm = CreateVector(n);
    XdbTp = CreateVector(n);
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

        Xdbi = OswinIsotherm(o, awi, T*(1-dev));
        setvalV(XdbTm, i, Xdbi);

        Xdbi = OswinIsotherm(o, awi, T*(1+dev));
        setvalV(XdbTp, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k0 = o->k0*(1-dev);
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbk0m, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k0 = o->k0*(1+dev);
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbk0p, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k1 = o->k1*(1-dev);
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbk1m, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k1 = o->k1*(1+dev);
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbk1p, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->n0 = o->n0*(1-dev);
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbn0m, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->n0 = o->n0*(1+dev);
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbn0p, i, Xdbi);
        DestroyOswinData(o);
        
        o = CreateOswinData();
        o->n1 = o->n1*(1-dev);
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbn1m, i, Xdbi);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->n1 = o->n1*(1+dev);
        Xdbi = OswinIsotherm(o, awi, T);
        setvalV(Xdbn1p, i, Xdbi);
        DestroyOswinData(o);
    }

    dev = dev * 100;
    hdr = (char*) calloc(sizeof(char), 500);
    sprintf(hdr, "Xdb,D,T-%.0f%%,T+%.0f%%,k0-%.0f%%,k0+%.0f%%,k1-%.0f%%,k1+%.0f%%,n0-%.0f%%,n0+%.0f%%,n1-%.0f%%,n1+%.0f%%\n",
            dev, dev, dev, dev, dev,
            dev, dev, dev, dev, dev);
    out = CatColVector(12, aw, Xdb, XdbTm, XdbTp, Xdbk0m, Xdbk0p, Xdbk1m, Xdbk1p, Xdbn0m, Xdbn0p, Xdbn1m, Xdbn1p);

    //mtxprntfilehdr(out, "output.csv", "aw,Xdb,T-10%,T+10%,k0-10%,k0+10%,k1-10%,k1+10%,n0-10%,n0+10%,n1-10%,n1+10%\n");
    mtxprntfilehdr(out, "output.csv", hdr);

    printf("Deviation\tL2 Norm\n");
    printf("T-%.0f%%:\t\t%g\n", dev, l2errornorm(aw, Xdb, XdbTm));
    printf("T+%.0f%%:\t\t%g\n", dev, l2errornorm(aw, Xdb, XdbTp));
    printf("k0-%.0f%%:\t\t%g\n", dev, l2errornorm(aw, Xdb, Xdbk0m));
    printf("k0+%.0f%%:\t\t%g\n", dev, l2errornorm(aw, Xdb, Xdbk0p));
    printf("k1-%.0f%%:\t\t%g\n", dev, l2errornorm(aw, Xdb, Xdbk1m));
    printf("k1+%.0f%%:\t\t%g\n", dev, l2errornorm(aw, Xdb, Xdbk1p));
    printf("n0-%.0f%%:\t\t%g\n", dev, l2errornorm(aw, Xdb, Xdbn0m));
    printf("n0+%.0f%%:\t\t%g\n", dev, l2errornorm(aw, Xdb, Xdbn0p));
    printf("n1-%.0f%%:\t\t%g\n", dev, l2errornorm(aw, Xdb, Xdbn1m));
    printf("n1+%.0f%%:\t\t%g\n", dev, l2errornorm(aw, Xdb, Xdbn1p));

    return;
}

