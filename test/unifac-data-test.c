#include "../unifac/unifac.h"

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
    UnifacAddMolec(s, m, .7);

    DestroyVector(ids);
    DestroyVector(count);

    ids = CreateVector(2);
    setvalV(ids, 0, 10);
    setvalV(ids, 1, 14);
    count = CreateVector(2);
    setvalV(count, 0, 1);
    setvalV(count, 1, 2);

    m = UnifacCreateMolec(ids, count, d);
    UnifacAddMolec(s, m, .3);

    UnifacPrintSolution(s);

    return 0;
}

