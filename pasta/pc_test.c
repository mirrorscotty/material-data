#include "pasta.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    vector *Xdb, *Pc, *aw;
    matrix *out;
    oswin *o;

    int L = 500,
        i;

    double Xmin = .0001, Xmax = .3, T = 273+60;

    Xdb = CreateVector(L);
    Pc = CreateVector(L);
    aw = CreateVector(L);
    o = CreateOswinData();

    for(i = 0; i < L; i++) {
        setvalV(Xdb, i, Xmin + i*((Xmax-Xmin)/L));
    }
    for(i = 0; i < L; i++) {
        setvalV(aw, i, OswinInverse(o, valV(Xdb, i), T));
        setvalV(Pc, i, pore_press(valV(Xdb, i), T));
    }

    out = CatColVector(3, Xdb, aw, Pc);

    mtxprntfile(out, "Pc.csv");

    return 0;
}

