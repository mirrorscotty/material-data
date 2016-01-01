#include "matrix.h"
#include "material-data.h"

#include <stdlib.h>

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

int main(int argc, char *argv[])
{
    vector *Xdb, *strain, *stress, *Pc, *phi, *aw;
    double Xmin = .005,
           Xmax = .33,
           Xi = .33,
           T,
           n = 100,
           J;
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
    stress = CreateVector(n);
    phi = CreateVector(n);
    aw = CreateVector(n);
    Pc = CreateVector(n);
    for(i=0; i<n; i++) {
        J = CreepLookupJ0("creep-333K.csv", T, valV(Xdb, i))
            + CreepLookupJ1("creep-333K.csv", T, valV(Xdb, i))
            + CreepLookupJ2("creep-333K.csv", T, valV(Xdb, i));
        setvalV(strain, i, linear_strain(Xi, valV(Xdb, i), T));
        setvalV(stress, i,
                valV(strain, i)/J);
        setvalV(Pc, i, EffPorePress(Xi, valV(Xdb, i), T));
        setvalV(aw, i, OswinInverse(o, valV(Xdb, i), T));
        setvalV(phi, i, solidfrac(Xi, T, valV(strain, i)));
    }

    output = CatColVector(6, Xdb, strain, stress, Pc, aw, phi);
    mtxprntfilehdr(output, "stress.csv", "Xdb,strain,stress,Pc,aw,solidfrac\n");
}

