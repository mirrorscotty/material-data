/**
 * @file diff-test.c
 * Testing file for diffusivity
 */

#include "isotherms.h"
#include "diffusivity.h"
#include "binding.h"
#include "matrix.h"
#include "choi-okos.h"
#include "constants.h"
#include "pasta.h"
#include <stdio.h>
#include <stdlib.h>

void diff_test() {
    printf("Diffusivity Test:\n");
    printf("T: %g Xdb: %g\t Deff = 82.8, %g\n", 105.+273, .21, DiffCh10(.21, 105+273));
    printf("T: %g Xdb: %g\t Deff = 35.8, %g\n", 71.+273, .21, DiffCh10(.21,71+273));
    printf("T: %g Xdb: %g\t Deff = 79.2, %g\n", 105.+273, .08, DiffCh10(.08, 105+273));
    printf("T: %g Xdb: %g\t Deff = 26.7, %g\n", 71.+273, .08, DiffCh10(.08, 71+273));
    printf("T: %g Xdb: %g\t Deff = 7.35, %g\n", 44.+273, .09, DiffCh10(.09,44+273));
    puts("");
    printf("T: %g Xdb: %g\t Deff = %g\n", 60+273., .04, DiffCh10(0.04,60+273));

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
    vector *X, *D40, *D80;
    //double D0 = 6.3910e-8;
    double D0 = 1;
    matrix *data;
    int i;

    X = linspaceV(.005, .18, 300);

    D40 = CreateVector(300);
    D80 = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(D40, i, DiffCh10GAB(valV(X, i), 40+273.15)/D0);
        setvalV(D80, i, DiffCh10GAB(valV(X, i), 80+273.15)/D0);
    }

    data = CatColVector(3, X, D40, D80);
    mtxprntfilehdr(data, "output.csv", "Xdb,D(40),D(80)\n");
}

void PlotEb()
{
    vector *X, *Eb30, *Eb40, *Eb50, *Eb60, *Eb70, *Eb80;
    matrix *data;
    int i;
    double A = 2.39e-4; /* Conversion to kcal from J */
    //oswin *d;
    gab *d;
    d = GABDATA();

    X = linspaceV(.01, .25, 300);

    Eb30 = CreateVector(300);
    Eb40 = CreateVector(300);
    Eb50 = CreateVector(300);
    Eb60 = CreateVector(300);
    Eb70 = CreateVector(300);
    Eb80 = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(Eb30, i, A*BindingEnergyGAB(d, valV(X, i), 30+273.15));
        setvalV(Eb40, i, A*BindingEnergyGAB(d, valV(X, i), 40+273.15));
        setvalV(Eb50, i, A*BindingEnergyGAB(d, valV(X, i), 50+273.15));
        setvalV(Eb60, i, A*BindingEnergyGAB(d, valV(X, i), 60+273.15));
        setvalV(Eb70, i, A*BindingEnergyGAB(d, valV(X, i), 70+273.15));
        setvalV(Eb80, i, A*BindingEnergyGAB(d, valV(X, i), 80+273.15));
    }
    data = CatColVector(7, X, Eb30, Eb40, Eb50, Eb60, Eb70, Eb80);
    mtxprntfilehdr(data, "EbGAB.csv", "Xdb,30C,40C,50C,60C,70C,80C\n");
}

void PlotEbLitchfield()
{
    vector *X, *Eb30, *Eb40, *Eb50, *Eb60, *Eb70, *Eb80;
    matrix *data;
    int i;
    double A = 2.39e-4; /* Conversion to kcal from J */
    //oswin *d;
    gab *d;
    d = GABDATA();

    X = linspaceV(.01, .25, 300);

    Eb30 = CreateVector(300);
    Eb40 = CreateVector(300);
    Eb50 = CreateVector(300);
    Eb60 = CreateVector(300);
    Eb70 = CreateVector(300);
    Eb80 = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(Eb30, i, A*BindingEnergyLitchfield(valV(X, i), 30+273.15));
        setvalV(Eb40, i, A*BindingEnergyLitchfield(valV(X, i), 40+273.15));
        setvalV(Eb50, i, A*BindingEnergyLitchfield(valV(X, i), 50+273.15));
        setvalV(Eb60, i, A*BindingEnergyLitchfield(valV(X, i), 60+273.15));
        setvalV(Eb70, i, A*BindingEnergyLitchfield(valV(X, i), 70+273.15));
        setvalV(Eb80, i, A*BindingEnergyLitchfield(valV(X, i), 80+273.15));
    }
    data = CatColVector(7, X, Eb30, Eb40, Eb50, Eb60, Eb70, Eb80);
    mtxprntfilehdr(data, "EbLitchfield.csv", "Xdb,30C,40C,50C,60C,70C,80C\n");
}

void PlotEbOswin()
{
    vector *X, *Eb30, *Eb40, *Eb50, *Eb60, *Eb70, *Eb80;
    matrix *data;
    int i;
    double A = 2.39e-4; /* Conversion to kcal from J */
    oswin *d;
    d = OSWINDATA();

    X = linspaceV(.01, .25, 300);

    Eb30 = CreateVector(300);
    Eb40 = CreateVector(300);
    Eb50 = CreateVector(300);
    Eb60 = CreateVector(300);
    Eb70 = CreateVector(300);
    Eb80 = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(Eb30, i, A*BindingEnergyOswin(d, valV(X, i), 30+273.15));
        setvalV(Eb40, i, A*BindingEnergyOswin(d, valV(X, i), 40+273.15));
        setvalV(Eb50, i, A*BindingEnergyOswin(d, valV(X, i), 50+273.15));
        setvalV(Eb60, i, A*BindingEnergyOswin(d, valV(X, i), 60+273.15));
        setvalV(Eb70, i, A*BindingEnergyOswin(d, valV(X, i), 70+273.15));
        setvalV(Eb80, i, A*BindingEnergyOswin(d, valV(X, i), 80+273.15));
    }
    data = CatColVector(7, X, Eb30, Eb40, Eb50, Eb60, Eb70, Eb80);
    mtxprntfilehdr(data, "EbOswin.csv", "Xdb,30C,40C,50C,60C,70C,80C\n");
}

void PlotEbHenderson()
{
    vector *X, *Eb30, *Eb40, *Eb50, *Eb60, *Eb70;
    matrix *data;
    int i;
    double A = 2.39e-4; /* Conversion to kcal from J */
    henderson *d;
    d = CreateHendersonData();

    X = linspaceV(.01, .25, 300);

    Eb30 = CreateVector(300);
    Eb40 = CreateVector(300);
    Eb50 = CreateVector(300);
    Eb60 = CreateVector(300);
    Eb70 = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(Eb30, i, A*BindingEnergyHenderson(d, valV(X, i), 30+273.15));
        setvalV(Eb40, i, A*BindingEnergyHenderson(d, valV(X, i), 40+273.15));
        setvalV(Eb50, i, A*BindingEnergyHenderson(d, valV(X, i), 50+273.15));
        setvalV(Eb60, i, A*BindingEnergyHenderson(d, valV(X, i), 60+273.15));
        setvalV(Eb70, i, A*BindingEnergyHenderson(d, valV(X, i), 70+273.15));
    }
    data = CatColVector(6, X, Eb30, Eb40, Eb50, Eb60, Eb70);
    mtxprntfile(data, "EbHenderson.csv");
}

void PlotXdb()
{
    vector *aw, *X40, *X55, *X71;
    matrix *data;
    int i;
    henderson *d;
    d = CreateHendersonData();

    aw = linspaceV(.005, .94, 300);

    X40 = CreateVector(300);
    X55 = CreateVector(300);
    X71 = CreateVector(300);

    for(i=0; i<len(aw); i++) {
        setvalV(X40, i, HendersonInverse(d, valV(aw, i), 40+273.15));
        setvalV(X55, i, HendersonInverse(d, valV(aw, i), 55+273.15));
        setvalV(X71, i, HendersonInverse(d, valV(aw, i), 71+273.15));
    }
    data = CatColVector(4, aw, X40, X55, X71);
    mtxprntfile(data, "XdbHenderson.csv");
}

void TestAw()
{
    henderson *d;

    vector *X, *Aw40, *Aw55, *Aw71;
    matrix *data;
    int i;

    d = CreateHendersonData();

    X = linspaceV(0.001, 0.9, 300);
    Aw40 = CreateVector(300);
    Aw55 = CreateVector(300);
    Aw71 = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(Aw40, i, HendersonIsotherm(d, valV(X, i), 40+273));
        setvalV(Aw55, i, HendersonIsotherm(d, valV(X, i), 55+273));
        setvalV(Aw71, i, HendersonIsotherm(d, valV(X, i), 71+273));
    }

    data = CatColVector(4, X, Aw40, Aw55, Aw71);
    mtxprntfile(data, "AwHenderson.csv");
}

void TestDCap()
{
    choi_okos *co;

    vector *X, *D401, *D551, *D711, *D402, *D552, *D712;
    matrix *data;
    int i;

    co = CreateChoiOkos(WATERCOMP);

    X = linspaceV(0.005, 0.5, 300);
    D401 = CreateVector(300);
    D551 = CreateVector(300);
    D711 = CreateVector(300);
    D402 = CreateVector(300);
    D552 = CreateVector(300);
    D712 = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(D401, i, CapillaryDiff(valV(X, i), 40+273));
        setvalV(D551, i, CapillaryDiff(valV(X, i), 55+273));
        setvalV(D711, i, CapillaryDiff(valV(X, i), 71+273));
        setvalV(D402, i, DiffCh10(valV(X, i), 40+273));
        setvalV(D552, i, DiffCh10(valV(X, i), 55+273));
        setvalV(D712, i, DiffCh10(valV(X, i), 71+273));
    }

    data = CatColVector(7, X, D401, D551, D711, D402, D552, D712);
    mtxprntfile(data, "Dcap.csv");
}

void CompareDiffXdb(double X)
{
    vector *T, *D10o, *D10g, *D10m, *Dz1, *Dz2, *Dhend;
    char *filename;
    matrix *out;
    int i;

    filename = (char*) calloc(sizeof(char), 80);

    T = linspaceV(273.15+25, 273.15+90, 300);
    D10o = CreateVector(300);
    D10g = CreateVector(300);
    D10m = CreateVector(300);
    Dz1 = CreateVector(300);
    Dz2 = CreateVector(300);
    Dhend = CreateVector(300);

    for(i=0; i<len(T); i++) {
        setvalV(D10o, i, DiffCh10(X, valV(T, i)));
        setvalV(D10g, i, DiffCh10GAB(X, valV(T, i)));
        setvalV(Dhend, i, DiffCh10Hend(X, valV(T, i)));
        setvalV(D10m, i, DiffCh10Mod(X, valV(T, i)));
        setvalV(Dz1, i, CapillaryDiff(X, valV(T, i)));
        setvalV(Dz2, i, CapDiff(X, valV(T, i)));
    }
    sprintf(filename, "Diffusivity%g.csv", X);
    out = CatColVector(7, T, Dhend, D10o, D10g, D10m, Dz1, Dz2);
    mtxprntfile(out, filename);
}

void CompareAllDiff(double T)
{
    vector *X, *D10o, *D10g, *D10m, *Dz1, *Dz2,
           *Dhend, *Dlitchfield, *Dwaananen,
           *Dachanta, *DachantaOrig, *DachantaFit;
    char *filename;
    matrix *out;
    int i;

    filename = (char*) calloc(sizeof(char), 80);
    T = T+273.15;

    X = linspaceV(0.005, 0.3, 300);
    D10o = CreateVector(300);
    D10g = CreateVector(300);
    D10m = CreateVector(300);
    Dz1 = CreateVector(300);
    Dz2 = CreateVector(300);
    Dhend = CreateVector(300);
    Dlitchfield = CreateVector(300);
    Dwaananen = CreateVector(300);
    Dachanta = CreateVector(300);
    DachantaOrig = CreateVector(300);
    DachantaFit = CreateVector(300);

    for(i=0; i<len(X); i++) {
        setvalV(D10o, i, DiffCh10(valV(X, i), T));
        setvalV(D10g, i, DiffCh10GAB(valV(X, i), T));
        setvalV(Dhend, i, DiffCh10Hend(valV(X, i), T));
        setvalV(D10m, i, DiffCh10Mod(valV(X, i), T));
        setvalV(Dz1, i, CapillaryDiff(valV(X, i), T));
        setvalV(Dz2, i, CapDiff(valV(X, i), T));
        setvalV(Dlitchfield, i, DiffLitchfield(valV(X, i), T));
        setvalV(Dwaananen, i, DiffWaananen(valV(X, i), T));
        setvalV(Dachanta, i, DiffAchanta(valV(X, i), T));
        setvalV(DachantaOrig, i, AchantaOrigDeff(valV(X, i), T));
        setvalV(DachantaFit, i, DiffAchantaFitted(valV(X, i), T));
    }
    sprintf(filename, "Diffusivity%gK.csv", T);
    out = CatColVector(12, X, Dhend, D10o, D10g, D10m, Dz1, Dz2,
            Dlitchfield, Dwaananen, Dachanta, DachantaOrig, DachantaFit);
    mtxprntfilehdr(out, filename, "Xdb,Henderson,Oswin,GAB,Modified,Zhu1,Zhu2,Litchfield,Waananen,Achanta,AchantaOrig,AchantaFit\n");
}

void CompareAllIsotherm(double T)
{
    vector *aw, *Singh, *Xiong, *XiongR, *Erbas, *Andrieu, *Bressani, *Waananen, *Griessman;
    vector *Che, *Kir;
    char *filename;
    matrix *out;
    int i;
    gab *dsingh, *derbas, *dandrieu, *dbressani, *dwaananen;
    oswin *dxiong, *dxiongr, *dgriessman;
    gab *dkir, *dche;

    dxiong = CreateOswinXiong();
    dxiongr = CreateOswinXiongR();
    dsingh = CreateGABSingh();
    derbas = CreateGABErbas();
    dandrieu = CreateGABAndrieu();
    dbressani = CreateGABData();
    dwaananen = CreateGABWaananen();
    dgriessman = CreateOswinGriessman();

    dkir = CreateGABPotatoKir();
    dche = CreateGABPotatoChemkhi();

    filename = (char*) calloc(sizeof(char), 80);
    T = T+273.15;

    aw = linspaceV(0.005, .99, 300);
    Xiong = CreateVector(300);
    XiongR = CreateVector(300);
    Singh = CreateVector(300);
    Erbas = CreateVector(300);
    Andrieu = CreateVector(300);
    Bressani = CreateVector(300);
    Waananen = CreateVector(300);
    Kir = CreateVector(300);
    Che = CreateVector(300);
    Griessman = CreateVector(300);

    for(i=0; i<len(aw); i++) {
        setvalV(Xiong, i, OswinIsotherm(dxiong, valV(aw, i), T));
        setvalV(XiongR, i, OswinIsotherm(dxiongr, valV(aw, i), T));
        setvalV(Singh, i, GABIsotherm(dsingh, valV(aw, i), T));
        setvalV(Erbas, i, GABIsotherm(derbas, valV(aw, i), T));
        setvalV(Andrieu, i, GABIsotherm(dandrieu, valV(aw, i), T));
        setvalV(Bressani, i, GABIsotherm(dbressani, valV(aw, i), T));
        setvalV(Waananen, i, GABIsotherm(dwaananen, valV(aw, i), T));
        setvalV(Kir, i, GABIsotherm(dkir, valV(aw, i), T));
        setvalV(Che, i, GABIsotherm(dche, valV(aw, i), T));
        setvalV(Griessman, i, OswinIsotherm(dgriessman, valV(aw, i), T));
    }
    sprintf(filename, "Isotherm%gK.csv", T);
    out = CatColVector(9, aw, Xiong,XiongR, Singh, Erbas, Andrieu, Bressani, Waananen, Griessman);
    mtxprntfilehdr(out, filename, "aw,Xiong,XiongR,Singh,Erbas,Andrieu,Bressani,Waananen,Griessman\n");
}

void Dgas(double P)
{
    vector *T, *D;
    matrix *out;
    int i;

    T = linspaceV(273.15, 373.15, 300);
    D = CreateVector(300);

    for(i=0; i<len(T); i++) {
        setvalV(D, i, VaporDiff(valV(T,i), P));
    }

    out = CatColVector(2, T, D);
    mtxprntfile(out, "GasDiff.csv");
}

matrix* DOswinVector(matrix *X, double T)
{
    double i;
    matrix *D;
    D = CreateMatrix(nRows(X), 1);
    for(i=0; i<nRows(X); i++)
        setval(D, DiffCh10(val(X, i, 0), T), i, 0);
    return D;
}

int main(int argc, char *argv[])
{
    matrix *data, *X, *D, *data1;
    double phi = POROSITY;
//    printf("Water Activity: %g\n", GABInverse(data, .5, 20));
    //diff_test();
/*    printf("D:\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n",
            CapillaryDiff(d, data, .1787, 60+273.15),
            CapillaryDiff(d, data, .1156, 60+273.15),
            CapillaryDiff(d, data, .0862, 60+273.15),
            CapillaryDiff(d, data, .069, 60+273.15),
            CapillaryDiff(d, data, .0577, 60+273.15),
            CapillaryDiff(d, data, .0489, 60+273.15),
            0.,
            CapillaryDiff(d, data, .0285, 60+273.15)); */

//    CompareAllDiff(105);
//    CompareAllDiff(71);
//    CompareAllDiff(55);
//    CompareAllDiff(44);
//
//      CompareDiffXdb(.05);
//      CompareDiffXdb(.1);
//      CompareDiffXdb(.2);
    CompareAllDiff(44);
    CompareAllDiff(55);
    CompareAllDiff(71);
    CompareAllDiff(105);
    CompareAllDiff(40);
    CompareAllDiff(70);
    CompareAllDiff(85);


    PlotEb();
    PlotEbLitchfield();
    PlotEbOswin();
    CompareAllIsotherm(40);
    CompareAllIsotherm(60);
    CompareAllIsotherm(80);
    PlotDeff();

/*    printf("44 = %g, 55 = %g, 71 = %g, 105 = %g\n",
            VaporDiff(44+273.15, 101325),
            VaporDiff(55+273.15, 101325),
            VaporDiff(71+273.15, 101325),
            VaporDiff(105+273.15, 101325)); */
    //Dgas(101325);

    /*printf("XsOswin = %g, XsSat = %g, phi = %g\n",
            OswinIsotherm(o, .95, 273.15+25),
            mdb_wat_sat(phi, 273.15+25),
            phi); */
    /* printf("Dfit @ 0C = %g, Dfit @ 50C = %g, Dfit @ 100C = %g\n",
            SelfDiffWater(273.15),
            SelfDiffWater(50+273.15),
            SelfDiffWater(373.15)); */

//    data = mtxloadcsv("kF.csv", 0);
//    X = ExtractColumn(data, 0);
//    D = DOswinVector(X, 60+273.15);
//    data1 = AugmentMatrix(data, D);
//    mtxprntfile(data1, "D-kF.csv");
    //diff_test();
    //PlotEbHenderson();
    //TestAw();
    //PlotXdb();
    return 0;
}

