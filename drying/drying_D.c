#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "drying_D.h"
#include "regress.h"

/* Used in Dv */
/* Used in Dl */
double Do, K, Ea;
/* Used in Cg */
double Cg0, EaCg;
/* Used in m0 */
double m00, Eam0;

double R, K0, EaK;
/*
int main(int argc, char *argv[])
{
    set_default_globals();

    * Test EVERYTHING *
    printf("Kk:  %f\n", Kk(280));
    printf("m0:  %f\n", m0(280));
    printf("Cg:  %f\n", Cg(280));
    printf("Aw:  %f\n", Aw(280, .5));
    printf("Ebf: %f\n", Ebf(280, 1));
    printf("Dl:  %e\n", Dl(280, .5));
    printf("Dv:  %e\n", Dv(280, .5));
    printf("D:   %e\n", drying_D(280, .5, .5));
    * Validates against Matlab Code *

    return 0;
}
*/

void set_default_globals()
{
    R = 8.314;
    K = 1032.558;
    
    Do = 2.607e-7;
    Ea = 25.31;
    
    Cg0 = 56899.584;
    EaCg = 22954.327;

    m00 = 0.0013897;
    Eam0 = -21951.174;

    K0 = 2.1554732;
    EaK = 2898.6332;
}

double drying_D(double T, double M, double por)
{
    double DV, DL;
    DV = Dv(T, M);
    DL = Dl(T, M);
    return por*DV + (1-por)*DL;
}

double Dv(double T, double M)
{
    double a, b, Eav, D0v;
    /* Global variables: R, K */
    matrix *Dveff, *Tvinv, *z;
    /* [0, 1/25, 1/42] */
    Tvinv = ParseMatrix("[0, 0.04, 2.38095e-2]");
    Dveff = ParseMatrix("[0.22e-4, 0.26e-4, 0.288e-4]");
    Map(Dveff, &log);

    z = polyfit(Tvinv, Dveff, 1);
    a = val(z, 1, 0);
    b = val(z, 0, 0);

    D0v = exp(b);
    Eav = -a*R;

    DestroyMatrix(Dveff);
    DestroyMatrix(Tvinv);
    DestroyMatrix(z);

    /* Needs work... lots of work */
    return D0v*exp(-Eav/(R*T))*(K*exp(-Ebf(T,M)/(R*T))/(1+K*exp(-Ebf(T,M)/(R*T))));
}

double Dl(double T, double M)
{
    double x;
    /* Global variables: Do, Ea, K, R */
    x = Do * exp(-Ea/(R*T)) * ( (K*exp(-Ebf(T,M)/(R*T))) / (1+K*exp(-Ebf(T,M)/(R*T))) );
    return x;
}

double Ebf(double T, double x)
{
    matrix *mc, *Eb, *lnEb, *mcc, *p;
    int i, j, Ebzeroindex, Ebmaxindex;
    int Mcfinal;
    double a, b, bb, ln, TT;
    /* Global variables: R */
    Mcfinal = 30;
    Eb = CreateMatrix(1, Mcfinal);

    mc = linspace(1, Mcfinal, Mcfinal);
    for(i=0; i<Mcfinal; i++) {
        ln = log(Aw(T+0.001, i+1)/Aw(T,i+1));
        TT = (1/T) - (1/(T+0.001));
        setval(Eb, R*ln/TT, 0, i);
        //printf("%f*%f/%f = %e\n", R, ln, TT, R*ln/TT);
    }
    Ebmaxindex = FindMaxIndex(Eb);
    Ebzeroindex = FindFirstZeroIndex(Eb);

    //mtxprnt(Eb);

    lnEb = CreateMatrix(1, Ebzeroindex-Ebmaxindex);
    mcc = CreateMatrix(1, Ebzeroindex-Ebmaxindex);
    if(!lnEb->rows) {
        fprintf(stderr, "Terminated.\n");
        exit(1);
    }
    if(!mcc->rows) {
        fprintf(stderr, "Terminated.\n");
        exit(1);
    }
    
    /* This might be backwards */
    for(i=Ebmaxindex, j=0; i<Ebzeroindex; i++, j++) {
        setval(lnEb, log(val(Eb, 0, i)), 0, j);
        setval(mcc, val(mc, 0, i), 0, j);
    }
    //mtxprnt(mcc);
    //mtxprnt(lnEb);
    p = polyfit(mcc, lnEb, 2);
    a = val(p, 2, 0);
    b = val(p, 1, 0);
    bb = exp(b);

    /* Free stuff! */
    DestroyMatrix(mc);
    DestroyMatrix(Eb);
    DestroyMatrix(lnEb);
    DestroyMatrix(mcc);
    DestroyMatrix(p);

    return bb*exp(x*a);
}

int FindFirstZeroIndex(matrix *A)
{
    int i;
    for(i=0; i<nCols(A); i++)
        if(val(A, 0, i) == 0)
            return i;
    return i;
}

int FindMaxIndex(matrix *A)
{
    double max;
    int maxindex, i;

    max = 0;

    for(i=0; i<nCols(A); i++) {
        if(val(A, 0, i) > max) {
            maxindex = i;
            max = val(A, 0, i);
        }
    }

    return maxindex;
}

double Aw(double T, double xw)
{
    double x, xx;
    xx = .5*(xw*Cg(T) - 2*xw - m0(T)*Cg(T) + pow(pow(xw,2)*pow(Cg(T),2)-2*xw*pow(Cg(T),2)*m0(T)+4*xw*m0(T)*Cg(T)+pow(m0(T),2)*pow(Cg(T),2), .5))/xw/(Cg(T)-1)/Kk(T);
    if(xx>1)
        x = 1;
    else
        x = xx;
    
    return x;
}

double Cg(double Tt)
{
    double x;
    /* Global varibles: R, Cg0, EaCg */
    x = Cg0*exp(-EaCg/(R*Tt));
    return x;
}

double m0(double Tt)
{
    double x;
    /* Global variables: R, m00, Eam0 */
    x = m00*exp(-Eam0/(R*Tt));
    return x;
}

double Kk(double Tt)
{
    double x;
    /* Global variables used: R, K0, EaK */
    x = K0*exp(-EaK/(R*Tt));
    return x;
}

