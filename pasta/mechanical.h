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

typedef struct {
    double J0;
    double *J;

    double *tau;

    double mu0;

    int n;
} burgers;

maxwell* CreateMaxwell();
void DestroyMaxwell(maxwell*);
double ReducedTime(maxwell*, double, double, double);
double MaxwellRelax(maxwell*, double, double, double);
double MaxwellCreep(maxwell*, double, double, double);
double MeanRelaxTime(maxwell*);
double MaxwellStress(maxwell*, double, double (*)(double), double (*)(double),
        double (*)(double));
double pore_press(double, double);

burgers* CreateBurgers();
burgers* CreateBurgers1();
void DestroyBurgers();
double BurgersCreep(burgers*, double, double, double);
double DBurgersCreep(burgers*, double, double, double);

#endif

