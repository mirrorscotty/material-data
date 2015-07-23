#include <math.h>
#include <complex.h>

double complex laplace(double complex (*f)(double complex, void*),
                       double complex s, void *param)
{
    const int DefaultIntegralN = 5000;
    double du = .5 / DefaultIntegralN,
           u = 0,
           limit = 1.0 - 1e-10;
    double complex y = -f(1e-10, param) / 2.0;

    while(u < limit) {
        u += du;
        y += 2.0 * cpow(u, s-1) * f(clog(u), param);
        u += du;
        y += cpow(u, s-1) * f(clog(u), param);
    }

    return 2.0 * y * du / 3.0;
}

