/**
 * @file diffusivity.h
 */

#ifndef DIFFUSIVITY_H
#define DIFFUSIVITY_H

#include "isotherms.h"

#define OSWINDATA CreateOswinXiong
#define GABDATA CreateGABAndrieu

double CapillaryDiff(double, double);
double CapDiff(double, double);
double CapillaryPressure(oswin*, double, double);
double DiffCh10(double, double);
double DiffCh10GAB(double, double);
double DiffCh10Mod(double, double);
double VaporDiff(double, double);
double VaporDiffCh10(double, double);
double SelfDiffWater(double);

#endif

