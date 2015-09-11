/**
 * @file salt-gsp.c
 * Water activities over saturated salt solutions according to Greenspan (1977)
 */

#include "material-data.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
    oswin *o;
    double T, Tc, awLiCl, awKAc, awMgCl2, awK2CO3, awNaCl, awKCl, awKNO3, awK2SO4, XdbTg;
    gordontaylor *gt;

    o = CreateOswinData();
    gt = GTSemolina();

    if(argc != 2) {
        puts("Usage: salt-calc <T>");
        puts("T: Temperature [C]");
        exit(0);
    }

    T = atof(argv[1]) + 273.15;
    Tc = T - 273.15;

    awLiCl = 11.2323 - 0.00824245 * Tc - 0.214890e-3 * Tc*Tc;
    awKAc = 22.4388 + 0.156288 * Tc + 0.612268e-2 * Tc*Tc;
    awMgCl2 = 33.6686 - 0.00797397 * Tc - 0.108988e-2 * Tc*Tc;
    awK2CO3 = 43.1315 - 0.00147523 * Tc;
    awNaCl = 75.5164 + 0.0398321 * Tc - 0.265459e-2 * Tc*Tc + 0.2848e-4 * Tc*Tc*Tc;
    awKCl = 88.6190 - 0.193340 * Tc + 0.899706e-3 * Tc*Tc;
    awKNO3 = 96.3361 + 0.0112371 * Tc - 0.484514e-2 * Tc*Tc;
    awK2SO4 = 98.7792 - 0.0590502 * Tc;
    XdbTg = GordonTaylorInv(gt, T);

    printf("T = %g\n", T);

    printf("Salt\tRH\t\tXeq\n");
    printf("LiCl\t%g%%\t%g\n", awLiCl, OswinIsotherm(o, awLiCl/100, T));
    printf("KAc\t%g%%\t%g\n", awKAc, OswinIsotherm(o, awKAc/100, T));
    printf("MgCl2\t%g%%\t%g\n", awMgCl2, OswinIsotherm(o, awMgCl2/100, T));
    printf("K2CO3\t%g%%\t%g\n", awK2CO3, OswinIsotherm(o, awK2CO3/100, T));
    printf("NaCl\t%g%%\t%g\n", awNaCl, OswinIsotherm(o, awNaCl/100, T));
    printf("KCl\t%g%%\t%g\n", awKCl, OswinIsotherm(o, awKCl/100, T));
    printf("KNO3\t%g%%\t%g\n", awKNO3, OswinIsotherm(o, awKNO3/100, T));
    printf("K2SO4\t%g%%\t%g\n", awK2SO4, OswinIsotherm(o, awK2SO4/100, T));

    printf("Tg\t%g%%\t%g\n", OswinInverse(o, XdbTg, T)*100, XdbTg);

    return 0;
}

