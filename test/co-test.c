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
        puts("Usage: co-test <T>");
        puts("T: Temperature [C]");
        exit(0);
    }
    T = atof(argv[1]);

    co = CreateChoiOkos(WATERCOMP);
    printf("rho_w = %g kg/m^3\n", rho(co, T));
    DestroyChoiOkos(co);
    co = CreateChoiOkos(PASTACOMP);
    printf("rho_s = %g\n kg/m^3", rho(co, T));
    DestroyChoiOkos(co);

    return 0;
}

