/**
 * @file diffusivity.h
 */

#ifndef DIFFUSIVITY_H
#define DIFFUSIVITY_H

#include "../isotherms/isotherms.h"
#include "matrix.h"

#define OSWINDATA CreateOswinXiong
#define GABDATA CreateGABAndrieu

typedef struct{
    double D0;
    double Ea;
    double K;
    double Ebf;
} DiffXiongData;

DiffXiongData* CreateDefaultXiongData();
void DestroyXiongData(DiffXiongData *d);

double CapillaryDiff(double, double);
double CapDiff(double, double);
double CapillaryPressure(oswin*, double, double);
double DiffCh10new(DiffXiongData*, oswin*, double, double);
double DiffCh10(double, double);
double DiffCh10GAB(double, double);
double DiffCh10Hend(double, double);
double DiffCh10Mod(double, double);
double VaporDiff(double, double);
double VaporDiffCh10(double, double);
double SelfDiffWater(double);

double DiffLitchfield(double, double);
double DiffWaananen(double, double);
double DiffChemkhi(double, double);

double DiffAchanta(double, double);
double AchantaOrigDeff(double, double);
double AchantaDiffModel(matrix*, matrix*);

#endif

