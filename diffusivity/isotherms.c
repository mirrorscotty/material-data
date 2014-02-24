#include <math.h>

#include "diff_data.h"
#include "isotherms.h"

double OswinIsotherm(oswin *dat, double aw, double T)
{
    double Xdb;
    Xdb = (dat->k0 + dat->k1*T) * pow(aw/(1-aw), (dat->n0 + dat->n1*T));
    return Xdb;
}

double OswinInverse(oswin *dat, double X, double T)
{
    double aw;
    T = T-273.15;
    aw = pow(X/(dat->k0 + dat->k1*T), 1/(dat->n0 + dat->n1*T))
        / (pow(X/(dat->k0 + dat->k1*T), 1/(dat->n0 + dat->n1*T)) + 1);
    return aw;
}

double OswinDawDx(oswin *dat, double X, double T)
{
    double D, A;
    A = pow(X/(dat->k0 + dat->k1*T), 1/(dat->n0 + dat->n1*T));
    D = A/(X*(dat->n0 + dat->n1*T)*pow(A+1, 2));
    return D;
}

double BindingEnergyOswin(oswin *dat, double X, double T)
{
    double Eb, h, R, dawd1T;
    h = .00001;
    R = 8.314;

    dawd1T = (OswinInverse(dat, X, T+h) - OswinInverse(dat, X, T-h))/(2*h) * (-T*T);

    Eb = -R*dawd1T;

    return Eb;
}

double GABIsotherm(gab *dat, double aw, double T)
{
    double Xdb, Xm, C, k;
    
    Xm = dat->m0*exp(dat->dHm/T);
    C = dat->C0*exp(dat->dHC/T);
    k = dat->k0*exp(dat->dHk/T);

    Xdb = Xm * (C*k*aw)/((1-k*aw)*(1-k*aw+C*k*aw));
    return Xdb;
}

double GABInverse(gab *dat, double X, double T)
{
    double xm, c, k, A, B, C, y, aw;
    T = T-273.15;

    xm = dat->m0*exp(dat->dHm/T);
    c = dat->m0*exp(dat->dHC/T);
    k = dat->k0*exp(dat->dHk/T);

    A = 1/(k*c);
    B = (k-2)/k;
    C = (c-k*c)/k;
    y = X/xm;

    aw = (-1*sqrt(-4*A*C*y*y + B*B*y*y - 2*B*y + 1) - B*y + 1)/(2*C*y);

    return aw;
}

double BindingEnergyGAB(gab *dat, double X, double T)
{
    double Eb, h, R, dawd1T;
    h = .00001;
    R = 8.314;

    dawd1T = (GABInverse(dat, X, T+h) - GABInverse(dat, X, T-h))/(2*h) * (-T*T);

    Eb = -R*dawd1T;

    return Eb;
}


