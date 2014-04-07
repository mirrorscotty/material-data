#ifndef ISOTHERMS_H
#define ISOTHERMS_H

/* Values for the Oswin isotherm correlation */
struct _oswin {
    double k0;
    double k1;
    double n0;
    double n1;
};
typedef struct _oswin oswin;

/* Values for the GAB isotherm equation */
struct _gab {
    double m0;
    double C0;
    double k0;
    double dHm;
    double dHk;
    double dHC;
};
typedef struct _gab gab;

oswin* CreateOswinData();
oswin* CreateOswinXiong();
void DestroyOswinData(oswin*);
gab* CreateGABData();
gab* CreateGABErbas();
gab* CreateGABAndrieu();

double OswinIsotherm(oswin*, double, double);
double OswinInverse(oswin*, double, double);
double OswinDawDx(oswin*, double, double);

double GABIsotherm(gab*, double, double);
double GABInverse(gab*, double, double);

double BindingEnergyGAB(gab*, double, double);
double BindingEnergyOswin(oswin*, double, double);

#endif

