#include "mechanical.h"
#include <stdlib.h>

mechdat* CreateMechDat(double T, double M, double P)
{
    mechdat *d;
    d = (mechdat*) calloc(sizeof(mechdat), 1);
    d->T = T;
    d->M = M;
    d->P = P;
    return d;
}

void DestroyMechDat(mechdat *d)
{
    free(d);
}

