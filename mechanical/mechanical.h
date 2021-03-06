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

typedef struct {
    double J0b;
    double *Jb;

    double J0e;
    double *Je;

    double *tau;

    double mu0b;
    double mu0e;

    int n;
} burgerse;

typedef struct {
    double T;
    double M;
    double P;
} mechdat;

maxwell* CreateMaxwell();
maxwell* CreateMaxwellZhu();
void DestroyMaxwell(maxwell*);
double TimeShift(maxwell*, double, double);
double ReducedTime(maxwell*, double, double, double);
double MaxwellRelax(maxwell*, double, double, double);
double DMaxwellRelax(maxwell*, double, double, double);
double MaxwellRelaxLaura(double, double, double);
double DMaxwellRelaxLaura(double, double, double);
double MaxwellCreep(maxwell*, double, double, double);
double DMaxwellCreep(maxwell*, double, double, double);
double MaxwellCreepLaura(double, double, double);
double DMaxwellCreepLaura(double, double, double);

double LCummingsCreep(double, double, double, double);
double DLCummingsCreep(double, double, double, double);
double LLauraCreep(double, double, double, double);
double DLLauraCreep(double, double, double, double);
double LGinaRelax(double, double, double, double);
double DLGinaRelax(double, double, double, double);

double MaxwellCreepConverted(double, double, double);
double DMaxwellCreepConverted(double, double, double);

double MeanRelaxTime(maxwell*);
double MaxwellStress(maxwell*, double, double (*)(double), double (*)(double),
        double (*)(double));
double pore_press(double, double);
double pore_press_exp(double, double);
double pore_press_gab(double, double);

burgers* CreateBurgers();
burgerse* CreateBurgersE();
burgers* CreateBurgers1();
void DestroyBurgers();
void DestroyBurgersE();
double BurgersCreep(burgers*, double, double, double);
double DBurgersCreep(burgers*, double, double, double);
double BurgersECreep(burgerse*, double, double, double, double);
double DBurgersECreep(burgerse*, double, double, double, double);

double poisson(double, double);
double porosity(double, double, double, double);
double solidfrac(double, double, double);

mechdat* CreateMechDat(double, double, double);
void DestroyMechDat(mechdat*);

double CreepLookupJ0(char*, double, double);
double CreepLookupJ1(char*, double, double);
double CreepLookupJ2(char*, double, double);
double CreepLookupTau1(char*, double, double);
double CreepLookupTau2(char*, double, double);

#endif

