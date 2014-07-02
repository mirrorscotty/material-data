#ifndef MECHANICAL_H
#define MECHANICAL_H

typedef struct {
    double *E;
    double *tau;

    double aT0;
    double T0;

    double aM0;
    double M0;

    int n;
} maxwell;

maxwell* CreateMaxwell();
void DestroyMaxwell(maxwell*);
double ReducedTime(maxwell*, double, double, double);
double MaxwellRelax(maxwell*, double, double, double);
double MaxwellCreep(maxwell*, double, double, double);
double MeanRelaxTime(maxwell*);
double MaxwellStress(maxwell*, double, double (*)(double), double (*)(double),
        double (*)(double));
double pore_press(double, double);

#endif

