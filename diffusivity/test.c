#include "diff_data.h"
#include "isotherms.h"
#include "diffusivity.h"
#include "matrix/matrix.h"
#include "../choi-okos/choi-okos.h"
#include <stdio.h>

void GAB_test()
{
    gab *d;
    double T, X;
    d = CreateGABData();
    printf("GAB Test\n");
    T = 293; X = .2;
    printf("1: %g\t2: %g\t3: %g\n", GABInverse(d, X, T), GABInverse2(d, X, T), GABInverse3(d, X, T));
    T = 303; X = .2;
    printf("1: %g\t2: %g\t3: %g\n", GABInverse(d, X, T), GABInverse2(d, X, T), GABInverse3(d, X, T));
    T = 293; X = .5;
    printf("1: %g\t2: %g\t3: %g\n", GABInverse(d, X, T), GABInverse2(d, X, T), GABInverse3(d, X, T));
    T = 343; X = .2;
    printf("1: %g\t2: %g\t3: %g\n", GABInverse(d, X, T), GABInverse2(d, X, T), GABInverse3(d, X, T));
    return;
}

void diff_test() {
    printf("Diffusivity Test:\n");
    printf("T: %g Xdb: %g\t Deff = 82.8, %g\n", 105.+273, .21, DiffCh10(.21, 105+273));
    printf("T: %g Xdb: %g\t Deff = 35.8, %g\n", 71.+273, .21, DiffCh10(.21,71+273));
    printf("T: %g Xdb: %g\t Deff = 79.2, %g\n", 105.+273, .08, DiffCh10(.08, 105+273));
    printf("T: %g Xdb: %g\t Deff = 26.7, %g\n", 71.+273, .08, DiffCh10(.08, 71+273));
    printf("T: %g Xdb: %g\t Deff = 7.35, %g\n", 44.+273, .09, DiffCh10(.09,44+273));
}

void Eb_test() {
    oswin *d;
    d = CreateOswinData();
    printf("Binding Energy Test:\n");
    printf("T: %g, Xdb: %g\t Eb = 15029, %g\n", 25+273.15, .1, BindingEnergyOswin(d, .1, 25+273.15));
    printf("T: %g, Xdb: %g\t Eb = 15029, %g\n", 50+273.15, .1, BindingEnergyOswin(d, .1, 50+273.15));

    printf("T: %g, Xdb: %g\t Eb = 5381, %g\n", 25+273.15, .15, BindingEnergyOswin(d, .15, 25+273.15));
    printf("T: %g, Xdb: %g\t Eb = 5381, %g\n", 50+273.15, .15, BindingEnergyOswin(d, .15, 50+273.15));

    printf("T: %g, Xdb: %g\t Eb = 1962, %g\n", 25+273.15, .2, BindingEnergyOswin(d, .2, 25+273.15));
    printf("T: %g, Xdb: %g\t Eb = 1962, %g\n", 50+273.15, .2, BindingEnergyOswin(d, .2, 50+273.15));
}

void PlotDeff()
{
    vector *X, *D40, *D55, *D71;
    double D0 = 6.3910e-8;
    FILE *fp;
    matrix *data;
    int i;

    X = linspaceV(.005, .5, 300);

    D40 = CreateVector(300);
    D55 = CreateVector(300);
    D71 = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(D40, i, DiffCh10(valV(X, i), 40+273.15)/D0);
        setvalV(D55, i, DiffCh10(valV(X, i), 55+273.15)/D0);
        setvalV(D71, i, DiffCh10(valV(X, i), 71+273.15)/D0);
    }

    fp = fopen("Deff.csv", "w");
    fprintf(fp, "\"Xdb\",\"D40\",\"D55\",\"D71\"\n");
    for(i=0; i<len(X); i++)
        //fprintf(fp, "%f,%g,%g,%g\n", valV(X,i), valV(D40,i)/D0, valV(D55,i)/D0, valV(D71,i)/D0);
        fprintf(fp, "%f,%g,%g,%g\n", valV(X,i), valV(D40,i), valV(D55,i), valV(D71,i));
    fclose(fp);
}

void PlotEb()
{
    vector *X, *Eb30, *Eb40, *Eb50, *Eb60, *Eb70;
    matrix *data;
    int i;
    oswin *d;
    d = CreateOswinXiong();

    X = linspaceV(.01, .25, 300);

    Eb30 = CreateVector(300);
    Eb40 = CreateVector(300);
    Eb50 = CreateVector(300);
    Eb60 = CreateVector(300);
    Eb70 = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(Eb30, i, BindingEnergyOswin(d, valV(X, i), 30+273.15));
        setvalV(Eb40, i, BindingEnergyOswin(d, valV(X, i), 40+273.15));
        setvalV(Eb50, i, BindingEnergyOswin(d, valV(X, i), 50+273.15));
        setvalV(Eb60, i, BindingEnergyOswin(d, valV(X, i), 60+273.15));
        setvalV(Eb70, i, BindingEnergyOswin(d, valV(X, i), 70+273.15));
    }
    data = CatColVector(6, X, Eb30, Eb40, Eb50, Eb60, Eb70);
    mtxprntfile(data, "EbGAB.csv");
}

void PlotXdb()
{
    vector *aw, *X40, *X55, *X71;
    matrix *data;
    int i;
    gab *d;
    d = CreateGABAndrieuK();

    aw = linspaceV(.005, .94, 300);

    X40 = CreateVector(300);
    X55 = CreateVector(300);
    X71 = CreateVector(300);

    for(i=0; i<len(aw); i++) {
        setvalV(X40, i, GABIsotherm(d, valV(aw, i), 25+273.15));
        setvalV(X55, i, GABIsotherm(d, valV(aw, i), 55+273.15));
        setvalV(X71, i, GABIsotherm(d, valV(aw, i), 90+273.15));
    }
    data = CatColVector(4, aw, X40, X55, X71);
    mtxprntfile(data, "XdbGAB.csv");
}

void TestDCap()
{
    diff_data *d;
    oswin *o;
    choi_okos *co;

    vector *X, *D401, *D551, *D711, *D402, *D552, *D712;
    matrix *data;
    int i;

    d = CreateDiffData();
    o = CreateOswinData();
    co = CreateChoiOkos(0, 0, 0, 0, 0, 1, 0);
    
    X = linspaceV(0.005, 0.5, 300);
    D401 = CreateVector(300);
    D551 = CreateVector(300);
    D711 = CreateVector(300);
    D402 = CreateVector(300);
    D552 = CreateVector(300);
    D712 = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(D401, i, CapillaryDiff(d, o, valV(X, i), 40+273));
        setvalV(D551, i, CapillaryDiff(d, o, valV(X, i), 55+273));
        setvalV(D711, i, CapillaryDiff(d, o, valV(X, i), 71+273));
        setvalV(D402, i, DiffCh10(valV(X, i), 40+273));
        setvalV(D552, i, DiffCh10(valV(X, i), 55+273));
        setvalV(D712, i, DiffCh10(valV(X, i), 71+273));
    }

    printf("Xs: 40-%g 55-%g 71-%g\n", OswinIsotherm(o, .95, 40+273), OswinIsotherm(o, .95, 55+273), OswinIsotherm(o, .95, 71+273));
    printf("rho: 40-%g 55-%g 71-%g\n", rho(co, 40+273), rho(co, 55+273), rho(co, 71+273));

    data = CatColVector(7, X, D401, D551, D711, D402, D552, D712);
    mtxprntfile(data, "Dcap.csv");
}

int main(int argc, char *argv[])
{
    oswin *data;
    diff_data *d;
    data = CreateOswinData();
    d = CreateDiffData();

//    printf("Water Activity: %g\n", GABInverse(data, .5, 20));
    //diff_test();
    PlotXdb();
//
/*    printf("D:\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n",
            CapillaryDiff(d, data, .1787, 60+273.15),
            CapillaryDiff(d, data, .1156, 60+273.15),
            CapillaryDiff(d, data, .0862, 60+273.15),
            CapillaryDiff(d, data, .069, 60+273.15),
            CapillaryDiff(d, data, .0577, 60+273.15),
            CapillaryDiff(d, data, .0489, 60+273.15),
            0.,
            CapillaryDiff(d, data, .0285, 60+273.15)); */

    return 0;
}

