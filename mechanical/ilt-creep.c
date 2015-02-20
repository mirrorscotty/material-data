#include "material-data.h"
#include "matrix.h"
#include "inv-laplace.h"
#include <complex.h>
#include <math.h>

double complex LMaxwellCreep(maxwell *m,
                             double complex s,
                             double T,
                             double M)
{
    double complex sum = 0;
    int i;
    double A = TimeShift (m, T, M);

    for(i=0; i<m->n; i++)
        sum += m->E[i]/(s+A/m->tau[i]);

    return sum;
}

double complex DLMaxwellCreep(maxwell *m,
                             double complex s,
                             double T,
                             double M)
{
    double complex sum = 0;
    int i;
    double A = TimeShift (m, T, M);

    for(i=0; i<m->n; i++)
        sum += -1*m->E[i]*A/m->tau[i] * 1/(s+A/m->tau[i]);

    return sum;
}

double complex LMaxwellLauraCreep(double complex s,
                                  double T,
                                  double M)
{
    double Ea, E1, E2, l1, l2, cg;

    Ea = 68.18*(1/(1+exp((M-250.92*exp(-0.0091*T))/2.19))+0.078);
    E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238));
    E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    l1 = 7;
    l2 = 110;
    Ea *= 1e6;
    E1 *= 1e6;
    E2 *= 1e6;
 
    return Ea/s + E1/(s+1/l1) + E2/(s+1/l2);
}

double complex DLMaxwellLauraCreep(double complex s,
                                   double T,
                                   double M)
{
    double E1, E2, l1, l2;

    E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238));
    E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    l1 = 7;
    l2 = 110;
    E1 *= 1e6;
    E2 *= 1e6;
 
    return -1*(E1/l1 * 1/(s+1/l1) + E2/l2 * 1/(s+1/l2));
}

double complex _LCummingsCreep(double complex s, void *params)
{
    maxwell *m;
    mechdat *d;
    double result, T, M;
           
    d = (mechdat*) params;
    m = CreateMaxwell();

    T = d->T;
    M = d->M;

    result =  1/(s*s*LMaxwellCreep(m, s, T, M));

    DestroyMaxwell(m);

    return result;
}

double complex _LLauraCreep(double complex s, void *params)
{
    mechdat *d;
    d = (mechdat*) params;
    double T = d->T,
           M = d->M;

    return 1/(s*s*LMaxwellLauraCreep(s, T, M));
}

double complex _DLCummingsCreep(double complex s, void *params)
{
    maxwell *m;
    mechdat *d;
    d = (mechdat*) params;
    double T = d->T,
           M = d->M;
    m = CreateMaxwell();

    return 1/(s*s*LMaxwellCreep(m, s, T, M));
}

double complex _DLLauraCreep(double complex s, void *params)
{
    mechdat *d;
    d = (mechdat*) params;
    double T = d->T,
           M = d->M;

    return 1/(s*s*LMaxwellLauraCreep(s, T, M));
}
double LCummingsCreep(double t, double T, double M, double P)
{
    vector *tv, *Gv;
    mechdat *param;
    double G;
    tv = CreateVector(1);
    setvalV(tv, 0, t);

    param = CreateMechDat(T, M, P);

    Gv = ilt_euler(&_LCummingsCreep, tv, 32, (void*)param);
    G = valV(Gv, 0);

    DestroyVector(Gv);
    DestroyVector(tv);
    DestroyMechDat(param);

    return G;
}
 
double LLauraCreep(double t, double T, double M, double P)
{
    vector *tv, *Gv;
    mechdat *param;
    double G;
    tv = CreateVector(1);
    setvalV(tv, 0, t);

    param = CreateMechDat(T, M, P);

    Gv = ilt_euler(&_LLauraCreep, tv, 32, (void*)param);
    G = valV(Gv, 0);

    DestroyVector(Gv);
    DestroyVector(tv);
    DestroyMechDat(param);

    return G;
}

double DLCummingsCreep(double t, double T, double M, double P)
{
    vector *tv, *Gv;
    mechdat *param;
    double G;
    tv = CreateVector(1);
    setvalV(tv, 0, t);

    param = CreateMechDat(T, M, P);

    Gv = ilt_euler(&_DLCummingsCreep, tv, 32, (void*)param);
    G = valV(Gv, 0);

    DestroyVector(Gv);
    DestroyVector(tv);
    DestroyMechDat(param);

    return G;
}
 
double DLLauraCreep(double t, double T, double M, double P)
{
    vector *tv, *Gv;
    mechdat *param;
    double G;
    tv = CreateVector(1);
    setvalV(tv, 0, t);

    param = CreateMechDat(T, M, P);

    Gv = ilt_euler(&_DLLauraCreep, tv, 32, (void*)param);
    G = valV(Gv, 0);

    DestroyVector(Gv);
    DestroyVector(tv);
    DestroyMechDat(param);

    return G;
}

