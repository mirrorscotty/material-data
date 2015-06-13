#include "material-data.h"
#include "matrix.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    vector *t, *Jcummings, *DJcummings, *Jlaura, *DJlaura, *Jlaura_ilt, *Jgina;
    matrix *out;
    burgerse *b;
    double ti, Jlaurai, Jginai, Jcummingsi, DJcummingsi, DJlaurai;
    int i, n=1000;
    double T=273+40, M=.1, P=.2e6;

    t = linspaceV(0, 1e8, n);
    Jgina = CreateVector(n);
    Jlaura = CreateVector(n);
    DJlaura = CreateVector(n);
    Jcummings = CreateVector(n);
    DJcummings = CreateVector(n);

    b = CreateBurgersE();

    for(i=0; i<n; i++) {
        ti = valV(t, i);
        Jginai = BurgersECreep(b, ti, T, M, P);
        Jlaurai = LLauraCreep(ti, T, M, P);
        DJlaurai = DLLauraCreep(ti, T, M, P);
        Jcummingsi = LCummingsCreep(ti, T, M, P);
        DJcummingsi = DLCummingsCreep(ti, T, M, P);

        setvalV(Jgina, i, Jginai);
        setvalV(Jlaura, i, Jlaurai);
        setvalV(DJlaura, i, DJlaurai);
        setvalV(Jcummings, i, Jcummingsi);
        setvalV(DJcummings, i, DJcummingsi);
    }

    out = CatColVector(4, t, Jcummings, Jlaura, Jgina);

    mtxprntfilehdr(out, "output.csv", "Time,Cummings,Rozzi,Bressani\n");

    DestroyBurgersE(b);
    DestroyVector(t);
    DestroyVector(Jcummings);
    DestroyVector(Jlaura);
    DestroyVector(Jgina);
    DestroyMatrix(out);
    return;
}


