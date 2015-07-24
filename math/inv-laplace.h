#ifndef INV_LAPLACE_H
#define INV_LAPLACE_H

#include <complex.h>

vector* ilt_euler(double complex (*)(double complex, void*), vector*, int, void*);
double complex laplace(double complex (*)(double complex, void*),
                       double complex, void *);
double l2errornorm(vector*, vector*, vector*);

#endif
