#ifndef INV_LAPLACE_H
#define INV_LAPLACE_H

#include <complex.h>

vector* ilt_euler(double complex (*)(double complex, void*), vector*, int, void*);

#endif