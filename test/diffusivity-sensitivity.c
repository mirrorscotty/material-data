#include "material-data.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../math/inv-laplace.h"

int main(int argc, char *argv[])
{
    vector *Xdb, *D, *DTm, *DTp, *Dk0m, *Dk0p, *Dk1m, *Dk1p, *Dn0m, *Dn0p, *Dn1m, *Dn1p, *DD0m, *DD0p, *DEam, *DEap, *DKm, *DKp, *DEbfm, *DEbfp;
    matrix *out;
    oswin *o;
    DiffXiongData *d;
    double Di, Xdbi, T, dev;
    int i, n=100;
    char *hdr;

    if(argc != 3) {
        printf("Usage:\n"
               "diffusivity-sensitivity <T> <dev>\n"
               "<T>\tTemperature (K)\n"
               "<dev>\tAmount to change each parameter by\n");
        exit(0);
    }

    T = atof(argv[1]);
    dev = atof(argv[2]);

    Xdb = linspaceV(0.005, .50, n);
    D = CreateVector(n);
    DTm = CreateVector(n);
    DTp = CreateVector(n);
    Dk0m = CreateVector(n);
    Dk0p = CreateVector(n);
    Dk1m = CreateVector(n);
    Dk1p = CreateVector(n);
    Dn0m = CreateVector(n);
    Dn0p = CreateVector(n);
    Dn1m = CreateVector(n);
    Dn1p = CreateVector(n);
    DD0m = CreateVector(n);
    DD0p = CreateVector(n);
    DEam = CreateVector(n);
    DEap = CreateVector(n);
    DKm = CreateVector(n);
    DKp = CreateVector(n);
    DEbfm = CreateVector(n);
    DEbfp = CreateVector(n);

    for(i=0; i<n; i++) {
        Xdbi = valV(Xdb, i);

        d = CreateDefaultXiongData();
        o = CreateOswinData();
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(D, i, Di);

        Di = DiffCh10new(d, o, Xdbi, T*(1-dev));
        setvalV(DTm, i, Di);

        Di = DiffCh10new(d, o, Xdbi, T*(1+dev));
        setvalV(DTp, i, Di);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k0 = o->k0*(1-dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(Dk0m, i, Di);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k0 = o->k0*(1+dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(Dk0p, i, Di);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k1 = o->k1*(1-dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(Dk1m, i, Di);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->k1 = o->k1*(1+dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(Dk1p, i, Di);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->n0 = o->n0*(1-dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(Dn0m, i, Di);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->n0 = o->n0*(1+dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(Dn0p, i, Di);
        DestroyOswinData(o);
        
        o = CreateOswinData();
        o->n1 = o->n1*(1-dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(Dn1m, i, Di);
        DestroyOswinData(o);

        o = CreateOswinData();
        o->n1 = o->n1*(1+dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(Dn1p, i, Di);
        DestroyOswinData(o);
        DestroyXiongData(d);

        o = CreateOswinData();

        d = CreateDefaultXiongData();
        d->D0 = d->D0*(1-dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(DD0m, i, Di);
        DestroyXiongData(d);

        d = CreateDefaultXiongData();
        d->D0 = d->D0*(1+dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(DD0p, i, Di);
        DestroyXiongData(d);

        d = CreateDefaultXiongData();
        d->Ea = d->Ea*(1-dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(DEam, i, Di);
        DestroyXiongData(d);

        d = CreateDefaultXiongData();
        d->Ea = d->Ea*(1+dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(DEap, i, Di);
        DestroyXiongData(d);

        d = CreateDefaultXiongData();
        d->K = d->K*(1-dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(DKm, i, Di);
        DestroyXiongData(d);

        d = CreateDefaultXiongData();
        d->K = d->K*(1+dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(DKp, i, Di);
        DestroyXiongData(d);

        d = CreateDefaultXiongData();
        d->Ebf = d->Ebf*(1-dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(DEbfm, i, Di);
        DestroyXiongData(d);

        d = CreateDefaultXiongData();
        d->Ebf = d->Ebf*(1+dev);
        Di = DiffCh10new(d, o, Xdbi, T);
        setvalV(DEbfp, i, Di);
        DestroyXiongData(d);
        DestroyOswinData(o);
    }

    out = CatColVector(20, Xdb, D, DTm, DTp, Dk0m, Dk0p, Dk1m, Dk1p, Dn0m, Dn0p, Dn1m, Dn1p, DD0m, DD0p, DEam, DEap, DKm, DKp, DEbfm, DEbfp);

    dev = dev * 100;
    hdr = (char*) calloc(sizeof(char), 500);
    sprintf(hdr, "Xdb,D,T-%.0f%%,T+%.0f%%,k0-%.0f%%,k0+%.0f%%,k1-%.0f%%,k1+%.0f%%,n0-%.0f%%,n0+%.0f%%,n1-%.0f%%,n1+%.0f%%,D0-%.0f%%,D0+%.0f%%,Ea-%.0f%%,Ea+%.0f%%,K-%.0f%%,K+%.0f%%,Ebf-%.0f%%,Ebf+%.0f%%\n",
            dev, dev, dev, dev, dev, dev,
            dev, dev, dev, dev, dev, dev,
            dev, dev, dev, dev, dev, dev);

    //mtxprntfilehdr(out, "output.csv", "Xdb,D,T-10%,T+10%,k0-10%,k0+10%,k1-10%,k1+10%,n0-10%,n0+10%,n1-10%,n1+10%,D0-10%,D0+10%,Ea-10%,Ea+10%,K-10%,K+10%,Ebf-10%,Ebf+10%\n");
    mtxprntfilehdr(out, "output.csv", hdr);

    printf("Input Change\tOutput Change\n");
    printf("T-%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, DTm));
    printf("T+%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, DTp));
    printf("k0-%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, Dk0m));
    printf("k0+%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, Dk0p));
    printf("k1-%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, Dk1m));
    printf("k1+%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, Dk1p));
    printf("n0-%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, Dn0m));
    printf("n0+%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, Dn0p));
    printf("n1-%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, Dn1m));
    printf("n1+%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, Dn1p));
    printf("D0-%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, DD0m));
    printf("D0+%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, DD0p));
    printf("Ea-%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, DEam));
    printf("Ea+%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, DEap));
    printf("K-%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, DKm));
    printf("K+%.0f%%:\t\t%g\n", dev, l2errornorm(Xdb, D, DKp));
    printf("Ebf-%.0f%%:\t%g\n", dev, l2errornorm(Xdb, D, DEbfm));
    printf("Ebf+%.0f%%:\t%g\n", dev, l2errornorm(Xdb, D, DEbfp));

    return;
}

