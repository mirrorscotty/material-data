/**
 * @file pc_test.c
 * Test file for capillary pressure.
 */
#include "pasta.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    vector *Xdb, *Pc1, *Pc2, *aw1, *aw2, *Eb1, *Eb2;
    matrix *out;
    oswin *o;

    int L = 500, /* Number of points to calculate */
        i;

    double Xmin = .0001, /* Minimum moisture content [kg/kg db] */
           Xmax = .5, /* Maximum moisture content [kg/kg db] */
           T1 = 273+40, /* Drying temperature [K] */
           T2 = 273+80; /* Drying temperature [K] */

    /* Create vectors */
    Xdb = CreateVector(L);
    Pc1 = CreateVector(L);
    Pc2 = CreateVector(L);
    aw1 = CreateVector(L);
    aw2 = CreateVector(L);
    Eb1 = CreateVector(L);
    Eb2 = CreateVector(L);
    /* Set up isotherm data */
    o = CreateOswinData();

    /* Fill in the moisture content vector */
    for(i = 0; i < L; i++)
        setvalV(Xdb, i, Xmin + i*((Xmax-Xmin)/L));
    for(i = 0; i < L; i++) {
        /* Water activity */
        setvalV(aw1, i, OswinInverse(o, valV(Xdb, i), T1));
        setvalV(aw2, i, OswinInverse(o, valV(Xdb, i), T2));
        /* Binding energy */
        setvalV(Eb1, i, BindingEnergyOswin(o, valV(Xdb, i), T1));
        setvalV(Eb2, i, BindingEnergyOswin(o, valV(Xdb, i), T2));
        /* Capillary pressure */
        setvalV(Pc1, i, pore_press(valV(Xdb, i), T1));
        setvalV(Pc2, i, pore_press(valV(Xdb, i), T2));
    }

    out = CatColVector(7, Xdb, aw1, aw2, Eb1, Eb2, Pc1, Pc2);

    mtxprntfilehdr(out, "output.csv", "Xdb,aw(40),aw(80),Eb(40),Eb(80),Pc(40),Pc(80)\n");

    return 0;
}

