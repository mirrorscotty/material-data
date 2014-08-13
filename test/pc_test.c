/**
 * @file pc_test.c
 * Test file for capillary pressure.
 */
#include "pasta.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    vector *Xdb, *Pc, *aw, *Eb;
    matrix *out;
    oswin *o;

    int L = 500, /* Number of points to calculate */
        i;

    double Xmin = .0001, /* Minimum moisture content [kg/kg db] */
           Xmax = .3, /* Maximum moisture content [kg/kg db] */
           T = 273+60; /* Drying temperature [K] */

    /* Create vectors */
    Xdb = CreateVector(L);
    Pc = CreateVector(L);
    aw = CreateVector(L);
    Eb = CreateVector(L);
    /* Set up isotherm data */
    o = CreateOswinData();

    /* Fill in the moisture content vector */
    for(i = 0; i < L; i++)
        setvalV(Xdb, i, Xmin + i*((Xmax-Xmin)/L));
    for(i = 0; i < L; i++) {
        /* Water activity */
        setvalV(aw, i, OswinInverse(o, valV(Xdb, i), T));
        /* Binding energy */
        setvalV(Eb, i, BindingEnergyOswin(o, valV(Xdb, i), T));
        /* Capillary pressure */
        setvalV(Pc, i, pore_press(valV(Xdb, i), T));
    }

    out = CatColVector(4, Xdb, aw, Eb, Pc);

    mtxprntfilehdr(out, "Pc.csv", "Xdb,aw,Eb,Pc\n");

    return 0;
}

