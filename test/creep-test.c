#include "material-data.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    vector *t, *Jgina, *Jlaura;
    matrix *out;
    burgerse *b;
    double ti, Jginai, Jlaurai;
    int i, n=1000;
    double T=298, M=.1, s=1e5;

    b = CreateBurgersE();

    t = linspaceV(0, 360, n);
    Jgina = CreateVector(n);
    Jlaura = CreateVector(n);

    for(i=0; i<n; i++) {
        ti = valV(t, i);
        Jginai = BurgersECreep(b, ti, T, M, s);
        Jlaurai = MaxwellCreepConverted(ti, T, M);

        setvalV(Jgina, i, Jginai);
        setvalV(Jlaura, i, Jlaurai);
    }

    out = CatColVector(3, t, Jgina, Jlaura);

    mtxprntfile(out, "output.csv");

    return;
}


