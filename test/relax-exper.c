#include "material-data.h"
#include "matrix.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double strain(double t)
{
    double str=.01, t0=10;
    if(t>t0)
        return str;
    else
        return 0;
}

int main(int argc, char *argv[])
{
    maxwell *m;
    vector *t, *Ecummings, *Elaura, *Egina;
    matrix *out;
    double ti, tj, Ecummingsi, Elaurai, Eginai, dstrain;
    int i, j, n=1000;
    double T,
           M,
           P=.2e6, /* Applied force (Pa) */
           tf=100, /* Final time */
           dt;
    char *outfile;

    if(argc < 3) {
        puts("Usage: relax-exper <T> <M>");
        puts("<T>\tTemperature\t\t(K)");
        puts("<M>\tMoisture Content\t(kg/kg db)");
        exit(0);
    }
    T = atof(argv[1]);
    M = atof(argv[2]);

    m = CreateMaxwell();

    t = linspaceV(0, tf, n);
    dt = valV(t, 1)-valV(t, 0);
    Ecummings = CreateVector(n);
    Elaura = CreateVector(n);
    Egina = CreateVector(n);

    for(i=0; i<n; i++) {
        ti = valV(t, i);

        Ecummingsi = 0;//MaxwellRelax(m, ti, T, M)*strain(ti);
        Elaurai = 0;//MaxwellRelaxLaura(ti, T, M)*strain(ti);
        Eginai = 0;//LGinaRelax(ti, T, M, P)*strain(ti);

        for(j=0; j<i; j++) {
            tj = valV(t, j);
            dstrain = (strain(tj+dt)-strain(tj))/dt;

            Ecummingsi += MaxwellRelax(m, ti, T, M)*dstrain*dt;
            Elaurai += MaxwellRelaxLaura(ti, T, M)*dstrain*dt;
            Eginai += LGinaRelax(ti, T, M, P)*dstrain*dt;
        }

        setvalV(Ecummings, i, Ecummingsi);
        setvalV(Elaura, i, Elaurai);
        setvalV(Egina, i, Eginai);
    }

    out = CatColVector(4, t, Ecummings, Elaura, Egina);

    outfile = (char*) calloc(30, sizeof(char));
    sprintf(outfile, "relax-%g-%g.csv", T, M);

    mtxprntfilehdr(out, outfile, "Time,Cummings,Rozzi,Bressani\n"); 

    DestroyMaxwell(m);
    DestroyVector(t);
    DestroyVector(Ecummings);
    DestroyVector(Elaura);
    DestroyVector(Egina);
    DestroyMatrix(out);
    return;
}

