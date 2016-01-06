#include "matrix.h"
#include "material-data.h"

#include <stdlib.h>
#include <math.h>

double EffPorePress(double Xi, double Xf, double T)
{
    double P = pore_press(Xf, T),
           P0 = pore_press(Xi, T);
    return -1*(P-P0);
}

double linear_strain(double Xi, double Xf, double T)
{
    double dX = Xf-Xi,
           b1 = .8839,
           b2 = 0.004264,
           b3 = -3.726e-5,
           b4 = 0.3704,
           b5 = 0.3149,
           b6 = -1.001e-5,
           LL0;
    T = T-273;

    LL0 = b1 + b2*T + b3*T*T + b4*dX + b5*dX*dX;
    return 1-LL0;
}

double stress_t_inf(double X, double T, double strain)
{
    double J;
    J = CreepLookupJ0("creep-333K.csv", T, X)
        + CreepLookupJ1("creep-333K.csv", T, X)
        + CreepLookupJ2("creep-333K.csv", T, X);

    return strain/J;
}

double stress_j0(double X, double T, double strain)
{
    double J;
    J = CreepLookupJ0("creep-333K.csv", T, X);

    return strain/J;
}

double stress_t(double X, double T, double strain, double t)
{
    double J, J0, J1, J2, tau1, tau2;
    J0 = CreepLookupJ0("creep-333K.csv", T, X);
    J1 = CreepLookupJ1("creep-333K.csv", T, X);
    J2 = CreepLookupJ2("creep-333K.csv", T, X);
    tau1 = CreepLookupTau1("creep-333K.csv", T, X);
    tau2 = CreepLookupTau2("creep-333K.csv", T, X);

    J = J0 + J1*(1-exp(-t/tau1)) + J2*(1-exp(-t/tau2)); 

    return strain/J;
}

int main(int argc, char *argv[])
{
    vector *Xdb, *strain, *Pc, *phi, *aw,
           *stressInf, *stressJ0, *stressT10, *stressT100,
           *stressT1000, *stressT10000;
    double Xmin = .075,
           Xmax = .35,
           Xi = .33,
           T,
           n = 100;
    matrix *output;
    int i;
    oswin *o;

    if(argc != 2) {
        puts("Usage: calc-stress <T>");
        puts("  <T>: Temperature in Kelvins");
        exit(0);
    }

    T = atof(argv[1]);

    o = CreateOswinData();

    Xdb = linspaceV(Xmin, Xmax, n);
    strain = CreateVector(n);
    phi = CreateVector(n);
    aw = CreateVector(n);
    Pc = CreateVector(n);

    stressInf = CreateVector(n);
    stressJ0 = CreateVector(n);
    stressT10 = CreateVector(n);
    stressT100 = CreateVector(n);
    stressT1000 = CreateVector(n);
    stressT10000 = CreateVector(n);

    for(i=0; i<n; i++) {
        setvalV(strain, i, linear_strain(Xi, valV(Xdb, i), T));
        setvalV(Pc, i, EffPorePress(Xi, valV(Xdb, i), T));
        setvalV(aw, i, OswinInverse(o, valV(Xdb, i), T));
        setvalV(phi, i, solidfrac(Xi, T, valV(strain, i)));

        setvalV(stressInf, i, stress_t_inf(valV(Xdb, i), T, valV(strain, i)));
        setvalV(stressJ0, i, stress_j0(valV(Xdb, i), T, valV(strain, i)));
        setvalV(stressT10, i, stress_t(valV(Xdb, i), T, valV(strain, i), 10));
        setvalV(stressT100, i, stress_t(valV(Xdb, i), T, valV(strain, i), 100));
        setvalV(stressT1000, i, stress_t(valV(Xdb, i), T, valV(strain, i), 1000));
        setvalV(stressT10000, i, stress_t(valV(Xdb, i), T, valV(strain, i), 10000));
    }

    output = CatColVector(11, Xdb, strain, Pc, aw, phi, stressInf, stressJ0,
                             stressT10, stressT100, stressT1000, stressT10000);
    mtxprntfilehdr(output, "stress.csv",
            "Xdb,strain,Pc,aw,solidfrac,tinf,j0,t10,t100,t1000,t10000\n");
}

