#ifndef ISOTHERMS_H
#define ISOTHERMS_H

#include "diff_data.h"

double OswinIsotherm(oswin*, double, double);
double OswinInverse(oswin*, double, double);
double OswinDawDx(oswin*, double, double);

double GABIsotherm(gab*, double, double);
double GABInverse(gab*, double, double);

double BindingEnergyGAB(gab*, double, double);
double BindingEnergyOswin(oswin*, double, double);

#endif

