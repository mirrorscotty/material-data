/**
 */

#include "material-data.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
    double T, rhos, rhow;
    choi_okos *co;

    if(argc != 2) {
        printf("Usage: %s <T>\n", argv[0]);
        puts("T: Temperature [K]");
        exit(0);
    }
    T = atof(argv[1]);

    co = CreateChoiOkos(WATERCOMP);
    puts("---- Water ----");
    printf("rho_w = %g kg/m^3\n", rho(co, T));
    printf("k_w = %g W/m-K\n", k(co, T));
    printf("Cp_w = %g J/kg-K\n", Cp(co, T));
    DestroyChoiOkos(co);
    co = CreateChoiOkos(PASTACOMP);
    puts("---- Pasta ----");
    printf("rho_s = %g kg/m^3\n", rho(co, T));
    printf("k_s = %g W/m-K\n", k(co, T));
    printf("Cp_s = %g J/kg-K\n", Cp(co, T));
    DestroyChoiOkos(co);
    co = CreateChoiOkos(GRAPEJUICECOMP);
    puts("---- Grape Juice ----");
    printf("rho = %g kg/m^3\n", rho(co, T));
    printf("k = %g W/m-K\n", k(co, T));
    printf("Cp = %g J/kg-K\n", Cp(co, T));
    DestroyChoiOkos(co);

    return 0;
}

