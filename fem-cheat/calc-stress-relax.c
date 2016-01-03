#include "matrix.h"
#include "material-data.h"
#include <math.h>
#include <stdlib.h>

double EaLaura(double M, double T)
{
    double Ea;
    M = M/(1+M); /* Convert from dry basis to wet basis */
    M *= 100; /* The moisture content should be converted to a percentage */

    Ea = 68.18*(1/(1+exp((M-250.92*exp(-0.0091*T))/2.19))+0.078);
    Ea *= 1e6;
    return Ea;
}

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
    vector *Xdb, *strain, *stress, *Pc, *phi, *aw, *xf;
    double Xmin = .05,
           Xmax = .29,
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
    stress = CreateVector(n);
    phi = CreateVector(n);
    aw = CreateVector(n);
    Pc = CreateVector(n);
    xf = CreateVector(n);
    for(i=0; i<n; i++) {
        setvalV(strain, i, linear_strain(Xi, valV(Xdb, i), T));
        setvalV(stress, i,
                valV(strain, i)*EaLaura(T, valV(Xdb, i)));
        setvalV(Pc, i, EffPorePress(Xi, valV(Xdb, i), T));
        setvalV(aw, i, OswinInverse(o, valV(Xdb, i), T));
        setvalV(xf, i, solidfrac(Xi, T, valV(strain, i)));
        setvalV(phi, i, porosity(Xi, valV(Xdb, i), T, valV(strain, i)));
    }

    output = CatColVector(7, Xdb, strain, stress, Pc, aw, xf, phi);
    mtxprntfilehdr(output, "stress.csv", "Xdb,strain,stress,Pc,aw,solidfrac,porosity\n");
}

