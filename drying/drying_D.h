#ifndef DRYING_D_H
#define DRYING_D_H
#include "regress.h"

void set_default_globals();
double drying_D(double, double, double);
double Dv(double, double);
double Dl(double, double);
double Ebf(double, double);
int FindFirstZeroIndex(matrix*);
int FindMaxIndex(matrix*);
double Aw(double, double);
double Cg(double);
double m0(double);
double Kk(double);

#endif
