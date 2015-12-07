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
oswin* CreateOswinXiongR();
oswin* CreateOswinGriessman();
void DestroyOswinData(oswin*);
gab* CreateGABData();
gab* CreateGABSingh();
gab* CreateGABErbas();
gab* CreateGABAndrieu();
gab* CreateGABWaananen();
henderson* CreateHendersonData();
void DestroyHendersonData(henderson*);

gab* CreateGABPotatoChemkhi();
gab* CreateGABPotatoKir();

double OswinIsotherm(oswin*, double, double);
double OswinInverse(oswin*, double, double);
double OswinDawDx(oswin*, double, double);
double OswinDlnawDx(oswin*, double, double);

double GABIsotherm(gab*, double, double);
double GABInverse(gab*, double, double);
double GABDawDx(gab*, double, double);
double GABDlnawDx(gab*, double, double);

double HendersonIsotherm(henderson*, double, double);
double HendersonInverse(henderson*, double, double);

#endif

