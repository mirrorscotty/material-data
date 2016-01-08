#include "matrix.h"
#include <stdio.h>
#include <math.h>

enum datacols {
    _M = 1,
    _J0 = 2,
    _J1 = 3,
    _J2 = 5,
    _TAU1 = 4,
    _TAU2 = 6
};

matrix *data;
int IsGo = 0;

double CreepLookup(char *file, double T, double M, int param)
{
    int i = 0;
    double M0, M1, y0, y1;
    //matrix *data;

    if(!IsGo) {
        data = mtxloadcsv(file, 1);
        IsGo = 1;
    }

    while(val(data, i, _M) < M)
        i++;

    M0 = val(data, i-1, _M);
    M1 = val(data, i, _M);
    y0 = val(data, i-1, param);
    y1 = val(data, i, param);

    //DestroyMatrix(data);
    return y0 + (y1 - y0) * (M-M0)/(M1-M0);
} 

double CreepLookupV2(char *file, double T, double M, int param)
{
    int i = 0;
    double M0, M1, y0, y1,
           Mmin, Mmax,
           n, i_e;

    if(!IsGo) {
        data = mtxloadcsv(file, 1);
        IsGo = 1;
    }

    n = nRows(data);
    Mmin = val(data, 0, _M);
    Mmax = val(data, n-1, _M);

    if(M < Mmin)
        M = 1.1*Mmin;
    else if(M > Mmax)
        M = 0.9*Mmax;

    i_e = (M - Mmin) * n/(Mmax-Mmin);

    M0 = val(data, floor(i_e), _M);
    M1 = val(data, ceil(i_e), _M);
    y0 = val(data, floor(i_e), param);
    y1 = val(data, ceil(i_e), param);

    return y0 + (y1 - y0) * (M-M0)/(M1-M0);
} 

double CreepLookupJ0(char *f, double T, double M)
{ return CreepLookupV2(f, T, M, _J0); }
double CreepLookupJ1(char *f, double T, double M)
{ return CreepLookupV2(f, T, M, _J1); }
double CreepLookupJ2(char *f, double T, double M)
{ return CreepLookupV2(f, T, M, _J2); }
double CreepLookupTau1(char *f, double T, double M)
{ return CreepLookupV2(f, T, M, _TAU1); }
double CreepLookupTau2(char *f, double T, double M)
{ return CreepLookupV2(f, T, M, _TAU2); }

