#ifndef DIFFUSIVITY_H
#define DIFFUSIVITY_H

#include "diff_data.h"

double CapillaryDiff(diff_data*, oswin*, double, double);
double CapDiff(diff_data*, oswin*, double, double);
double CapillaryPressure(diff_data*, oswin*, double, double);
double DiffCh10(double, double);
double DiffCh10GAB(double, double);

#endif

