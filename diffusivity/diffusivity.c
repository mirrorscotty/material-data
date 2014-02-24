#include "diff_data.h"
#include "isotherms.h"
#include <math.h>

//double CapillaryDiff(diff_data *data, oswin *o, double X, double T)
//{
//    double Dcap, DPcDSw;
//    choi_okos *co;
//    co = CreateChoiOkos(0, 0, 0, 0, 0, 1, 0);
//    aw = OswinInverse(o, X, T);
//    DawDX = OswinDawDx(o, X, T);
//    rhow = rho(co,T);
//    /* Estimate Xs by pretending that the material is saturated at aw=.95 */
//    Xs = OswinIsotherm(o, .95, T);
//    
//    DPcDSw = rho(co,T)*co->R*T/co->MW_wat * 1/aw * Xs * DawDX;
//
//    Dcap = -data->kw/(data->muw*data->phi) * DPcDSw;
//
//    DestroyChoiOkos(co);
//
//    return Dcap;
//}

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

