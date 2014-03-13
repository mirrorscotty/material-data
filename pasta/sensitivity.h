#ifndef SENSITIVITY_H
#define SENSITIVITY_H

char* fivevarcw(const char*, double (*)(double, double, double, double, double), double, double, double, double, double, double);
char* fivevarwv(const char*, double (*)(double, double, double, double, double), double, double, double, double, double, double);
char* fivevarphi(const char*, double (*)(double, double, double, double, double), double, double, double, double, double, double);
char* fivevarT(const char*, double (*)(double, double, double, double, double), double, double, double, double, double, double);
char* fivevarP(const char*, double (*)(double, double, double, double, double), double, double, double, double, double, double);

char* threevarcw(const char*, double (*)(double, double, double), double, double, double, double);
char* threevarphi(const char*, double (*)(double, double, double), double, double, double, double);
char* threevarT(const char*, double (*)(double, double, double), double, double, double, double);

char* onevar(const char*, double (*)(double), double, double);

#endif

