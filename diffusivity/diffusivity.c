#include "diff_data.h"
#include "../choi-okos/choi-okos.h"
#include "isotherms.h"
#include <math.h>

double CapillaryDiff(diff_data *data, oswin *o, double X, double T)
{
    double Dcap, DPcDSw, DawDX, aw, Xs, kw, Sw, Sr, fphi;
    choi_okos *co;

    /* Water activity from the Oswin isotherm model */
    aw = OswinInverse(o, X, T);

    /* Derivative of water activity with respect to d.b. moisture content */
    DawDX = OswinDawDx(o, X, T);

    /* Estimate Xs by pretending that the material is saturated at aw=.95 */
    Xs = OswinIsotherm(o, .95, T);

    /* Porosity factor */
    fphi = 1;

    /* Water permeability */
    Sr = .05/Xs; /* This value is completely made up */
    Sw = X/Xs;
    if(Sw>Sr)
        kw = data->kw*pow((Sw-Sr)/(1-Sr), 3)*fphi;
    else
        kw = 0;
    
    /* Derivative of capillary pressure with respect to water saturation.
     * Equation from Miranda and Silva 2005 */
    co = CreateChoiOkos(0, 0, 0, 0, 0, 1, 0);
    DPcDSw = -rho(co,T)*co->R*T/co->MW_wat * 1/aw * Xs * DawDX;
    DestroyChoiOkos(co);

    /* Capillary diffusivity from p104 of Zhu thesis */
    Dcap = -kw/(data->muw*data->phi) * DPcDSw;

    return Dcap; 
}

double CapDiff(diff_data *d, oswin *o, double X, double T)
{
    double D, DawDe, aw, vl, e, Xs, Sw, Sr, kw, fphi;
    choi_okos *co;

    /* Calculate the molar volume of water */
    co = CreateChoiOkos(0, 0, 0, 0, 0, 1, 0);
    vl = co->MW_wat/rho(co, T);
    DestroyChoiOkos(co);

    /* Estimate Xs by pretending that the material is saturated at aw=.95 */
    Xs = OswinIsotherm(o, .95, T);

    /* Porosity factor */
    fphi = 1;

    /* Water permeability */
    Sr = .02;
    Sw = X/Xs;
    if(Sw>Sr)
        kw = d->kw*pow((Sw-Sr)/(1-Sr), 3)*fphi;
    else
        kw = 0;
 
    /* Volume fraction of water */
    e = d->phi*X/Xs;

    /* Derivative of water activity with respect to volume fraction water */
    DawDe = OswinDawDx(o, X, T) * Xs/d->phi;
    
    /* Water activity from the Oswin isotherm model */
    aw = OswinInverse(o, X, T);

    D = e*kw/d->muw * d->R*T/vl * (log(aw) + e/aw * DawDe);
    return D;
}

double CapillaryPressure(diff_data *d, oswin *o, double X, double T)
{
    double aw, Pc;
    choi_okos *co;
    co = CreateChoiOkos(0, 0, 0, 0, 0, 1, 0);
    aw = OswinInverse(o, X, T);

    Pc = -rho(co,T)*co->R*T/co->MW_wat * log(aw);
    DestroyChoiOkos(co);
    return Pc;
}

double DiffCh10(double X, double T)
{
    oswin *dat;
    dat = CreateOswinData();

    double Deff;
    /* Source: Xiong et al (1991) */
    double D0 = 6.3910e-8;

    /* Source: Litchfield and Okos (1986) */
    double Ea = 25900;

    double K = 1032.6;
    double Eb = BindingEnergyOswin(dat, X, T);
    double R = 8.314;

    Deff = D0 * exp(-Ea/(R*T))
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

double DiffCh10GAB(double X, double T)
{
    gab *dat;
    dat = CreateGABData();

    double Deff;
    /* Source: Xiong et al (1991) */
    double D0 = 6.3910e-8;

    /* Source: Litchfield and Okos (1986) */
    double Ea = 25900;

    double K = 1032.6;
    double Eb = BindingEnergyGAB(dat, X, T);
    double R = 8.314;

    Deff = D0 * exp(-Ea/(R*T))
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

