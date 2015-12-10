#include "matrix.h"
#include "material-data.h"

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
    vector *Xdb, *strain, *stress, *Pc;
    double Xmin = .05,
           Xmax = .29,
           Xi = .33,
           T = 333,
           n = 100;
    matrix *output;
    int i;

    Xdb = linspaceV(Xmin, Xmax, n);
    strain = CreateVector(n);
    stress = CreateVector(n);
    Pc = CreateVector(n);
    for(i=0; i<n; i++) {
        setvalV(strain, i, linear_strain(Xi, valV(Xdb, i), T));
        setvalV(stress, i,
                valV(strain, i)/CreepLookupJ0("output.csv",
                                              T, valV(Xdb, i)));
        setvalV(Pc, i, EffPorePress(Xi, valV(Xdb, i), T));
    }

    output = CatColVector(4, Xdb, strain, stress, Pc);
    mtxprntfilehdr(output, "stress.csv", "Xdb,strain,stress,Pc\n");
}

