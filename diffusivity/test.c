#include "diff_data.h"
#include "isotherms.h"
#include "diffusivity.h"
#include "matrix/matrix.h"
#include <stdio.h>

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

    X = linspaceV(.05, .2, 100);
    D40 = CreateVector(100);
    D55 = CreateVector(100);
    D71 = CreateVector(100);

    for(i=0; i<len(X); i++) {
        setvalV(D40, i, DiffCh10(valV(X, i), 40+273.15));
        setvalV(D55, i, DiffCh10(valV(X, i), 55+273.15));
        setvalV(D71, i, DiffCh10(valV(X, i), 71+273.15));
    }

    fp = fopen("Deff.csv", "w");
    fprintf(fp, "\"Xdb\",\"D40\",\"D55\",\"D71\"\n");
    for(i=0; i<len(X); i++)
        fprintf(fp, "%f,%g,%g,%g\n", valV(X,i), valV(D40,i)/D0, valV(D55,i)/D0, valV(D71,i)/D0);
    fclose(fp);
}
    
    

int main(int argc, char *argv[])
{
    oswin *data;
    data = CreateOswinData();


//    printf("Water Activity: %g\n", GABInverse(data, .5, 20));
    //diff_test();
//    diff_test();
    PlotDeff();

    return 0;
}

