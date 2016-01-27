#include "matrix.h"
#include "material-data.h"
#include <math.h>

#define CINIT .4

double eff_pore_press_exp(double X, double T)
{
    double dX = X-CINIT,
           b0 = 1.011,
           b1 = 3.727e-4,
           b2 = 0.2407,
           //b3 = -1.080e-5,
           LL0, J, strain;
    T = T-273;

    LL0 = b0 + b1*T + b2*dX;
    strain = 1-LL0;

    if(T==333-273)
        J = CreepLookupJ0("creep-333K.csv", T, X)
            + CreepLookupJ1("creep-333K.csv", T, X)
            + CreepLookupJ2("creep-333K.csv", T, X);
    else if(T==313-273)
        J = CreepLookupJ0("creep-313K.csv", T, X)
            + CreepLookupJ1("creep-313K.csv", T, X)
            + CreepLookupJ2("creep-313K.csv", T, X);
    else if(T==353-273)
        J = CreepLookupJ0("creep-353K.csv", T, X)
            + CreepLookupJ1("creep-353K.csv", T, X)
            + CreepLookupJ2("creep-353K.csv", T, X);
    else
        J = 1;


    return -1*strain/J;
}

double eff_pore_press_exp_quadratic(double X, double T)
{
    double dX = X-CINIT,
           b1 = .8839,
           b2 = 0.004264,
           b3 = -3.726e-5,
           b4 = 0.3704,
           b5 = 0.3149,
           //b6 = -1.001e-5,
           LL0, J, strain;
    T = T-273;

    LL0 = b1 + b2*T + b3*T*T + b4*dX + b5*dX*dX;
    strain = 1-LL0;

    if(T==333-273)
        J = CreepLookupJ0("creep-333K.csv", T, X)
            + CreepLookupJ1("creep-333K.csv", T, X)
            + CreepLookupJ2("creep-333K.csv", T, X);
    else if(T==313-273)
        J = CreepLookupJ0("creep-313K.csv", T, X)
            + CreepLookupJ1("creep-313K.csv", T, X)
            + CreepLookupJ2("creep-313K.csv", T, X);
    else if(T==353-273)
        J = CreepLookupJ0("creep-353K.csv", T, X)
            + CreepLookupJ1("creep-353K.csv", T, X)
            + CreepLookupJ2("creep-353K.csv", T, X);
    else
        J = 1;


    return -1*strain/J;
}

double aw_from_pc(double X, double T)
{
    double Vm = 1.802e-5, /* m^3/mol */
           R = GASCONST,
           Pc = eff_pore_press_exp(X, T);
    return exp(Pc*Vm/(R*T) * 1/0.0612);
}

int main(int argc, char *argv[])
{
    matrix *output;
    vector *Xdb, *aw1, *aw2, *aw3,
           *awo1, *awo2, *awo3;
    double T1 = 40+273,
           T2 = 60+273,
           T3 = 80+273;
    int n = 100, i;
    oswin *o;

    Xdb = linspaceV(0.0005, .3, n);
    aw1 = CreateVector(n);
    aw2 = CreateVector(n);
    aw3 = CreateVector(n);
    awo1 = CreateVector(n);
    awo2 = CreateVector(n);
    awo3 = CreateVector(n);

    o = CreateOswinAndrieu();

    for(i=0; i<n; i++) {
        setvalV(aw1, i, aw_from_pc(valV(Xdb, i), T1));
        setvalV(aw2, i, aw_from_pc(valV(Xdb, i), T2));
        setvalV(aw3, i, aw_from_pc(valV(Xdb, i), T3));
        setvalV(awo1, i, OswinInverse(o, valV(Xdb, i), T1));
        setvalV(awo2, i, OswinInverse(o, valV(Xdb, i), T2));
        setvalV(awo3, i, OswinInverse(o, valV(Xdb, i), T3));
    }

    output = CatColVector(7, Xdb, aw1, aw2, aw3, awo1, awo2, awo3);
    mtxprntfilehdr(output, "output.csv", "Xdb,40C,60C,80C,Oswin 40C,Oswin 60C,Oswin 80C\n");
    return 0;
}

