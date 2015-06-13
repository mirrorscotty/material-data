#include "material-data.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    oswin *o;
    double T, Xdb;

    o = CreateOswinData();

    if(argc != 3) {
        puts("Usage: aw-calc <aw> <T>");
        puts("Xdb: Moisture content [decimal]");
        puts("T: Temperature [C]");
        exit(0);
    }

    Xdb = atof(argv[1]);
    T = atof(argv[2]) + 273.15;

    printf("Xdb = %g, T = %g\n", Xdb, T);

    printf("aw = %g\n", OswinInverse(o, Xdb, T));

    return 0;
}

