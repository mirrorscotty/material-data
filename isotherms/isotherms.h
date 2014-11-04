/**
 * @file isotherms.h
 */

#ifndef ISOTHERMS_H
#define ISOTHERMS_H

/**
 * Values for the Oswin isotherm correlation.
 */
typedef struct {
    double k0;
    double k1;
    double n0;
    double n1;
} oswin;

/**
 * Values for the GAB isotherm equation.
 */
typedef struct {
    double m0;
    double C0;
    double k0;
    double dHm;
    double dHk;
    double dHC;
} gab;

/**
 * Values for the modified Henderson isotherm from Litchfield 1992
 */
typedef struct {
    double A;
    double B;
    double C;
    double D;
} henderson;

oswin* CreateOswinData();
oswin* CreateOswinXiong();
void DestroyOswinData(oswin*);
gab* CreateGABData();
gab* CreateGABErbas();
gab* CreateGABAndrieu();
henderson* CreateHendersonData();
void DestroyHendersonData(henderson*);

double OswinIsotherm(oswin*, double, double);
double OswinInverse(oswin*, double, double);
double OswinDawDx(oswin*, double, double);

double GABIsotherm(gab*, double, double);
double GABInverse(gab*, double, double);

double HendersonIsotherm(henderson*, double, double);
double HendersonInverse(henderson*, double, double);

double BindingEnergyGAB(gab*, double, double);
double BindingEnergyOswin(oswin*, double, double);
double BindingEnergyHenderson(henderson*, double, double);

#endif

