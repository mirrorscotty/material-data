#include "material-data.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    vector *t, *Gcummings, *Glaura;
    matrix *out;
    maxwell *m;
    double ti, Gcummingsi, Glaurai;
    int i, n=1000;
    double T=298, M=.1;

    m = CreateMaxwell();

    t = linspaceV(0, 36e8, n);
    Gcummings = CreateVector(n);
    Glaura = CreateVector(n);

    for(i=0; i<n; i++) {
        ti = valV(t, i);
        Gcummingsi = MaxwellRelax(m, ti, T, M);
        Glaurai = MaxwellRelaxLaura(ti, T, M);

        setvalV(Gcummings, i, Gcummingsi);
        setvalV(Glaura, i, Glaurai);
    }

    out = CatColVector(3, t, Gcummings, Glaura);

    mtxprntfile(out, "output.csv");

    return;
}


