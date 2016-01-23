#include "material-data.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>

double exp_strain(double dX, double T)
{
    double b0 = .8839,
           b1 = .004264,
           b2 = -3.726e-5,
           b3 = .3704,
           b4 = .3149,
           b5 = -1.001e-5;
    T = T-273;
    return b0 + b1*T + b2*T*T + b3*dX + b4*dX*dX;
}

double exp_strain_lin(double dX, double T)
{
    double b0 = 1.011,
           b1 = -3.737e-4,
           b2 = 0.2407,
           b3 = -1.080e-5;
    T = T-273;
    return b0 + b1*T + b2*dX;
}

int main(int argc, char *argv[])
{
    vector *Xe, *phi1, *phi2, *phi3;
    matrix *out;
    int i, n=100;
    double Xei, phii, X0 = .33,
           T1=313, T2=333, T3=353,
           type=1;

    Xe = linspaceV(0, .3, n);
    phi1 = CreateVector(n);
    phi2 = CreateVector(n);
    phi3 = CreateVector(n);

    for(i=0; i<n; i++) {
        Xei = valV(Xe, i);

        phii = porosity(X0, Xei, T1, -1*(1-pow(exp_strain_lin(Xei-X0, T1),type)));
        //phii = exp_strain(Xei-X0, T1);
        setvalV(phi1, i, phii);
        phii = porosity(X0, Xei, T2, -1*(1-pow(exp_strain_lin(Xei-X0, T2),type)));
        //phii = exp_strain(Xei-X0, T2);
        setvalV(phi2, i, phii);
        phii = porosity(X0, Xei, T3, -1*(1-pow(exp_strain_lin(Xei-X0, T3),type)));
        //phii = exp_strain(Xei-X0, T3);
        setvalV(phi3, i, phii);

    }

    out = CatColVector(4, Xe, phi1, phi2, phi3);

    mtxprntfilehdr(out, "output.csv", "Xe,phi313,phi333,phi353\n");

    return;
}

