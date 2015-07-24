#include "material-data.h"
#include "matrix.h"
#include "../math/inv-laplace.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double complex GsNorm(double complex s, void *params)
{
    mechdat *d;
    d = (mechdat*) params;
    double T = d->T,
           M = d->M,
           Ea, E1, E2, l1, l2;
    double complex Gs;

    M = M/(1+M); /* Convert from dry basis to wet basis */
    M *= 100; /* The moisture content should be converted to a percentage */

    Ea = 68.18*(1/(1+exp((M-250.92*exp(-0.0091*T))/2.19))+0.078);
    E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238));
    E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    l1 = 7;
    l2 = 110;
    Ea *= 1e6;
    E1 *= 1e6;
    E2 *= 1e6;

    Gs = Ea/s + E1/(s+1/l1) + E2/(s+1/l2);

    return Gs;
}

double complex _LauraJs(double complex s, void *params)
{
    mechdat *d;
    d = (mechdat*) params;
    double T = d->T,
           M = d->M,
           Ea, E1, E2, l1, l2;
    double complex Gs;

    M = M/(1+M); /* Convert from dry basis to wet basis */
    M *= 100; /* The moisture content should be converted to a percentage */

    Ea = 68.18*(1/(1+exp((M-250.92*exp(-0.0091*T))/2.19))+0.078);
    E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238));
    E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    l1 = 7;
    l2 = 110;
    Ea *= 1e6;
    E1 *= 1e6;
    E2 *= 1e6;

    Gs = Ea/s + E1/(s+1/l1) + E2/(s+1/l2);

    return 1/(s*s*Gs);
}

double complex _Jt(double complex t, void *param)
{
    vector *tv, *Gv;
    double G;
    tv = CreateVector(1);
    setvalV(tv, 0, t);

    Gv = ilt_euler(&_LauraJs, tv, 32, param);
    G = valV(Gv, 0);

    DestroyVector(Gv);
    DestroyVector(tv);

    return G;
}

double complex _Gs(double complex s, void *param)
{
    double Js = laplace(&_Jt, s, param);
    return 1/(s*s*Js);
}

double Gt(double t, double T, double M, double P)
{
    vector *tv, *Gv;
    mechdat *param;
    double G;
    tv = CreateVector(1);
    setvalV(tv, 0, t);

    param = CreateMechDat(T, M, P);

    Gv = ilt_euler(&_Gs, tv, 32, (void*)param);
    G = valV(Gv, 0);

    DestroyVector(Gv);
    DestroyVector(tv);
    DestroyMechDat(param);

    return G;
}

int main(int argc, char *argv[])
{
    int n = 500,
        i;
    double T=333, M=.15;
    vector *Gnorm, *Gilt, *t;
    matrix *out;
    mechdat *param;

    param = CreateMechDat(T, M, 0);

    t = linspaceV(1, 1000, n);
    Gnorm = CreateVector(n);
    Gilt = CreateVector(n);

    for(i=0; i<n; i++) {
        printf("Point %d of %d\r", i, n);
        fflush(stdout);
        //setvalV(Gnorm, i, MaxwellRelaxLaura(valV(t, i), T, M));
        //setvalV(Gilt, i, Gt(valV(t, i), T, M, 0));
        setvalV(Gnorm, i, GsNorm(valV(t, i), param));
        setvalV(Gilt, i, _Gs(valV(t, i), param));
    }
    printf("Error: %g\n", l2errornorm(t, Gnorm, Gilt));

    out = CatColVector(3, t, Gnorm, Gilt);
    mtxprntfilehdr(out, "output.csv", "Time,G(normal),G(transformed)\n");

    return;
}

