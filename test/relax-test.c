#include "material-data.h"
#include "matrix.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    vector *t, *Gcummings, *Glaura, *Ggina, *DGgina;
    matrix *out;
    maxwell *m;
    double ti, Gcummingsi, Glaurai, Gginai, DGginai;
    int i, n=1000;
    double T=298, M=.3, P=.2e6;

    m = CreateMaxwell();

    t = linspaceV(0, 1e8, n);
    Gcummings = CreateVector(n);
    Glaura = CreateVector(n);
    Ggina = CreateVector(n);
    DGgina = CreateVector(n);

    for(i=0; i<n; i++) {
        ti = valV(t, i);
        Gcummingsi = MaxwellRelax(m, ti, T, M);
        Glaurai = MaxwellRelaxLaura(ti, T, M);
        Gginai = LGinaRelax(ti, T, M, P);
        DGginai = DLGinaRelax(ti, T, M, P);

        setvalV(Gcummings, i, Gcummingsi);
        setvalV(Glaura, i, Glaurai);
        setvalV(Ggina, i, Gginai);
        setvalV(DGgina, i, DGginai);
    }

    out = CatColVector(4, t, Gcummings, Glaura, Ggina);

    mtxprntfilehdr(out, "output.csv", "Time,Cummings,Rozzi,Bressani\n"); 

    DestroyMaxwell(m);
    DestroyVector(t);
    DestroyVector(Gcummings);
    DestroyVector(Glaura);
    DestroyVector(Ggina);
    DestroyMatrix(out);
    return;
}

