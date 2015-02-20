#include "material-data.h"
#include "matrix.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    vector *t, *Jcummings, *Jlaura, *Jlaura_ilt, *Jgina;
    matrix *out;
    burgerse *b;
    double ti, Jlaurai, Jginai, Jcummingsi;
    int i, n=1000;
    double T=298, M=.1, P=200000;

    t = linspaceV(0, 1e8, n);
    Jgina = CreateVector(n);
    Jlaura = CreateVector(n);
    Jcummings = CreateVector(n);

    b = CreateBurgersE();

    for(i=0; i<n; i++) {
        ti = valV(t, i);
        Jginai = BurgersECreep(b, ti, T, M, P);
        Jlaurai = LLauraCreep(ti, T, M, P);
        Jcummingsi = LCummingsCreep(ti, T, M, P);

        setvalV(Jgina, i, Jginai);
        setvalV(Jlaura, i, Jlaurai);
        setvalV(Jcummings, i, Jcummingsi);
    }

    //Jlaura_ilt = ilt_euler(&LLauraCreep, t, 32);

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


