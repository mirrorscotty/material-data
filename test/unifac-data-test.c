#include "../unifac/unifac.h"
#include <stdio.h>
#include <math.h>

double water_vapor_pressure(double T) {
    return exp(20.386 - 5132/T);
}

int main(int argc, char *argv[])
{
    vector *ids;
    vector *count;

    unifac_data *d;
    unifac_molec *m;
    unifac_solution *s;

    d = UnifacCreateData();
    UnifacLoadData(d, "unifac1.csv");
    UnifacLoadInteractions(d, "unifac2.csv");

    ids = CreateVector(1);
    setvalV(ids, 0, 6);
    count = CreateVector(1);
    setvalV(count, 0, 1);

    m = UnifacCreateMolec(ids, count, d);

    s = UnifacCreateSolution(d);
    UnifacAddMolec(s, m, 1);

    DestroyVector(ids);
    DestroyVector(count);

    ids = CreateVector(2);
    setvalV(ids, 0, 10);
    setvalV(ids, 1, 14);
    count = CreateVector(2);
    setvalV(count, 0, 1);
    setvalV(count, 1, 1);

    m = UnifacCreateMolec(ids, count, d);
    UnifacAddMolec(s, m, 0);

    UnifacPrintSolution(s);
    printf("gamma = %g\n", exp(ln_gamma(0, s, 298)));
    printf("pvap = %g\n", exp(ln_gamma(0, s, 298)) * water_vapor_pressure(298) * s->xi[0]);

    return 0;
}

