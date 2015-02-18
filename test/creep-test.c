#include "material-data.h"
#include "matrix.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

vector* ilt_euler(double complex (*)(double complex), vector*, int);

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

double complex LCummingsCreep(double complex s)
{
    maxwell *m;
    double T = 273+40,
           M = .08;
    m = CreateMaxwell();

    return 1/(s*s*LMaxwellCreep(m, s, T, M));
}

double complex LLauraCreep(double complex s)
{
    double T = 273+40,
           M = .08;

    return 1/(s*s*LMaxwellLauraCreep(s, T, M));
}

int main(int argc, char *argv[])
{
    vector *t, *Jcummings, *Jlaura, *Jlaura_ilt, *Jgina;
    matrix *out;
    burgerse *b;
    double ti, Jlaurai, Jginai;
    int i, n=1000;
    double T=298, M=.1, P=200000;

    t = linspaceV(0, 1e8, n);
    Jgina = CreateVector(n);
    Jlaura = CreateVector(n);

    b = CreateBurgersE();

    for(i=0; i<n; i++) {
        ti = valV(t, i);
        Jginai = BurgersECreep(b, ti, T, M, P);
        Jlaurai = MaxwellCreepLaura(ti, T, M);

        setvalV(Jgina, i, Jginai);
        setvalV(Jlaura, i, Jlaurai);
    }

    Jcummings = ilt_euler(&LCummingsCreep, t, 32);
    Jlaura_ilt = ilt_euler(&LLauraCreep, t, 32);

    out = CatColVector(5, t, Jcummings, Jlaura, Jlaura_ilt, Jgina);

    mtxprntfilehdr(out, "output.csv", "Time,Cummings,Rozzi,Rozzi ILT,Bressani\n"); 
    return;
}


