#include "material-data.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    oswin *o;
    double T, aw;

    o = CreateOswinData();

    if(argc != 3) {
        puts("Usage: aw-calc <aw> <T>");
        puts("aw: Water activity [decimal]");
        puts("T: Temperature [C]");
        exit(0);
    }

    aw = atof(argv[1]);
    T = atof(argv[2]) + 273.15;

    printf("aw = %g, T = %g\n", aw, T);

    printf("Xdb = %g\n", OswinIsotherm(o, aw, T));

    return 0;
}

