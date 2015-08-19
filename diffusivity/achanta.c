/**
 * @file achanta.c
 * All equations are taken from Achanta (1993)
 */

#include "material-data.h"
#include <math.h>

double viscosity(double T)
{
    double Tb = 273 + 100; /* Boiling point of water */
    /* Eq. 3.49 */
    return 7.93e-6 * exp(3.8*Tb/T);
}

double vapor_pressure(double T)
{
    /* Eq. 3.42 */
    return 1/760 * pow(10, 7.967 - 1668.0/(T-45.15));
}

double rho_v_w(oswin *d, double X, double T)
{
    double Rw = 461.52; /* (J/kg-K Gas Constant)/(Molar mass of water) */
    /* Eq. 3.41 */
    return vapor_pressure(T)*OswinInverse(d, X, T)/(Rw*T);
}

double rho_v(oswin *d, double X, double T, double P0)
{
    /* Eq. 3.40 */
    return 353.4*P0/T - 0.611*rho_v_w(d, X, T);
}

double DrhovwDx(oswin *d, double X, double T, double P)
{
    double Pv = vapor_pressure(T), /* Vapor pressure of water */
           aw = OswinInverse(d, X, T), /* Water activity */
           DawDx = OswinDawDx(d, X, T), /* d/dx[a_w] */
           Rw = 461.52, /* (J/kg-K Gas Constant)/(Molar mass of water) */
           result;

    /* d/dx[ <rhow>/<rhov> ], solved using Maxima */
    result = (353400000 * Pv * DawDx * P * Rw)
        /(124891560000*P*P*Rw*Rw - 431854800*Pv*aw*P*Rw + 373321*Pv*Pv*aw*aw);
    return result;
}

double permeability(oswin *d, double X, double T)
{
    double eta_w = viscosity(T), /* water viscosity */
           rho_w, /* water density */
           R = 8.314, /* Gas constant (kg/mol-K) */
           D0 = 1.78e-5, /* m^2/s */
           Ea = 8.81 * 4184, /* kcal/mol to J/mol */
           Dlfree = D0*exp(-Ea/(R*T)), /* Diffusivity of free water */
           Eb = BindingEnergyOswin(d, X, T), /* Binding energy */
           n = 2.1; /* Empirical parameter for normal pasta */
    choi_okos *co;

    co = CreateChoiOkos(WATERCOMP);
    rho_w = rho(co, T);
    DestroyChoiOkos(co);

    /* Eq. 3.53 */
    return Dlfree*exp(-Eb/(n*R*T)) * eta_w/(rho_w*R*T);
}

/**
 * Achanta's model for effective diffusivity. The diffusivity is a function of
 * temperature, moisture, and porosity, and was derived using a modified form
 * of Darcy's law.
 * @param X Moisture content [kg/kg db]
 * @param T Temperature [K]
 * @returns Effective diffusivity [m^2/s]
 */ 
double DiffAchanta(double X, double T)
{
    double P = 1, /* atm */
           K, /* permeability */
           rho_w, /* water density */
           rho_s, /* solid density */
           rhoV, /* <rho^v>^v */
           rhoVW, /* <rho^v_w>^v */
           DrhovwDX, /* d/dx [ <rho^v_w>^v/<rho^v>^v ] */
           epsilon = 0.07, /* Porosity(?) */
           DlnawDX, /* d/dx [ln a_w] */
           eta_w = viscosity(T), /* Water viscosity */
           R = 8.314, /* Gas constant kg/mol-K */
           Dvap = VaporDiff(T, P), /* Vapor diffusivity from BSL */
           Deff;
    oswin *d;

    d = CreateOswinData();

    choi_okos *co;
    co = CreateChoiOkos(PASTACOMP);
    rho_s = rho(co, T);
    DestroyChoiOkos(co);
    co = CreateChoiOkos(WATERCOMP);
    rho_w = rho(co, T);
    DestroyChoiOkos(co);

    //K = K0 * exp(-BindingEnergyOswin(d, X, T)/R*T),
    K = permeability(d, X, T);
    rhoV = rho_v(d,X,T,P),
    rhoVW = rho_v_w(d,X,T),
    DrhovwDX = DrhovwDx(d, X, T, P);


    /* Eq. 3.39 */
    Deff = K*rho_w/(eta_w*(1-epsilon)*rho_s) * (rho_w*R*T) * DlnawDX
        + rhoV*Dvap*DrhovwDX/(rho_w*(1+rhoVW/rhoV));

    DestroyOswinData(d);

    return Deff;
}

