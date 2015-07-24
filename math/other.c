#include "matrix.h"
#include <math.h>

double l2errornorm(vector *x, vector *f1, vector *f2)
{
    int i;
    double norm = 0, base = 0,
           dx, fa, fb;
    for(i=0; i<len(x)-1; i++) {
        dx = valV(x, i+1) - valV(x, i);
        fa = pow(valV(f1, i) - valV(f2, i), 2);
        fb = pow(valV(f1, i+1) - valV(f2, i+1), 2);
        norm += dx*(fb+fa);

        fa = pow(valV(f1, i), 2);
        fb = pow(valV(f1, i+1), 2);
        base += dx*(fb+fa);
    }
    return sqrt((.5*norm)/(.5*base));
}

