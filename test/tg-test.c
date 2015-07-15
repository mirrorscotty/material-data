#include "material-data.h"
#include "matrix.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    vector *Xdb, *aw, *Tg;
    matrix *out;
    oswin *o;
    gordontaylor *gt;
    double Tgi, awi, Xdbi;
    int i, n=100;

    o = CreateOswinData();
    gt = GTSemolina();

    Xdb = linspaceV(0, .5, n);
    aw = CreateVector(n);
    Tg = CreateVector(n);

    for(i=0; i<n; i++) {
        Xdbi = valV(Xdb, i);
        Tgi = GordonTaylor(gt, Xdbi);
        awi = OswinInverse(o, Xdbi, Tgi);

        setvalV(aw, i, awi);
        setvalV(Tg, i, Tgi);
    }

    out = CatColVector(3, aw, Xdb, Tg);

    mtxprntfilehdr(out, "output.csv", "aw,Xdb,Tg\n");
    printf("X40 = %g\nX60 = %g\nX80 = %g\n", GordonTaylorInv(gt, 40+273), GordonTaylorInv(gt, 60+273), GordonTaylorInv(gt, 80+273));

    return;
}


