/**
 * @file salt-calc.c
 * Outputs relative humidities over saturated salt solutions based on the
 * regressions equations from Rahman (2009).
 */

#include "material-data.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
    oswin *o;
    double T, awLiCl, awKAc, awMgCl2, awK2CO3, awNaCl, awKCl, XdbTg;
    gordontaylor *gt;

    o = CreateOswinData();
    gt = GTSemolina();

    if(argc != 2) {
        puts("Usage: salt-calc <T>");
        puts("T: Temperature [C]");
        exit(0);
    }

    T = atof(argv[1]) + 273.15;

    awLiCl = exp(500.95/T - 3.85);
    awKAc = exp(861.39/T - 4.33);
    awMgCl2 = exp(303.35/T - 2.13);
    awK2CO3 = exp(145/T - 1.3);
    awNaCl = exp(228.92/T - 1.04);
    awKCl = exp(367.58/T - 1.39);
    XdbTg = GordonTaylorInv(gt, T);

    printf("T = %g\n", T);

    printf("Salt\tRH\t\tXeq\n");
    printf("LiCl\t%g\t%g\n", awLiCl, OswinIsotherm(o, awLiCl, T));
    printf("KAc\t%g\t%g\n", awKAc, OswinIsotherm(o, awKAc, T));
    printf("MgCl2\t%g\t%g\n", awMgCl2, OswinIsotherm(o, awMgCl2, T));
    printf("K2CO3\t%g\t%g\n", awK2CO3, OswinIsotherm(o, awK2CO3, T));
    printf("NaCl\t%g\t%g\n", awNaCl, OswinIsotherm(o, awNaCl, T));
    printf("KCl\t%g\t%g\n", awKCl, OswinIsotherm(o, awKCl, T));

    printf("Tg\t%g\t%g\n", OswinInverse(o, XdbTg, T), XdbTg);

    return 0;
}

