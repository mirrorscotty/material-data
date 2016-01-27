#include "material-data.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    int i, n=100;
    double Xi, Di,
           T1=313, T2=333, T3=353,
           P = 101300,
           phi = .0612;
    vector *D1, *D2, *D3, *X;
    matrix *out;

    X = linspaceV(0.005, 0.3, n);
    D1 = CreateVector(n);
    D2 = CreateVector(n);
    D3 = CreateVector(n);

    for(i=0; i<n; i++) {
        Di = DeffModelTest(valV(X, i), T1, P, phi);
        setvalV(D1, i, Di);

        Di = DeffModelTest(valV(X, i), T2, P, phi);
        setvalV(D2, i, Di);

        Di = DeffModelTest(valV(X, i), T3, P, phi);
        setvalV(D3, i, Di);
    }

    out = CatColVector(4, X, D1, D2, D3);
    mtxprntfilehdr(out, "diffusivity.csv", "X,D313,D333,D353\n");

    return 0;
}

