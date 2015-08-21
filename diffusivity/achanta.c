/**
 * @file achanta.c
 * All equations are taken from Achanta (1993)
 */

#include "material-data.h"
#include "constants.h"
#include <math.h>
#include <stdio.h>

//#define R_CONST 8.3145 /* J/(kg-mol K) */

/** 
 * Viscosity as a function of temperature
 * @param T Temperature [K]
 * @returns Viscosity [Pa s]
 */
double viscosity(double T)
{
    double Tb = 273.15 + 100, /* Boiling point of water */
           result;

    /* Eq. 3.49 */
    result = 7.93e-6 * exp(3.8*Tb/T);

    return result;
}

/**
 * Vapor pressure of water as a function of temperature
 * @param T Temperature [K]
 * @returns Vapor Pressure [Pa]
 */
double vapor_pressure(double T)
{
    double result;
    /* Eq. 3.42 */
    result = 1./760. * pow(10, 7.967 - 1668.0/(T-45.15));
    result *= 101325; /* Convert from atm to Pa */
    
    return result;
}

/**
 * Volume-averaged mass density of the water in the vapor phase
 * \f[
 * <\rho^v_w>^v = \frac{p^0_v a_w}{R_w T}
 * \f]
 * @param d Oswin isotherm parameters
 * @param X Moisture content [kg/kg db]
 * @param T Temperature [K]
 * @returns Density [kg/m^3]
 */
double rho_v_w(oswin *d, double X, double T)
{
    double Rw = 461.52, /* (J/kg-K Gas Constant)/(Molar mass of water) */
           result; 
    /* Eq. 3.41 */
    result = vapor_pressure(T)*OswinInverse(d, X, T)/(Rw*T);

    return result;
}

/**
 * Density of the gas phase.
 * \f[
 * <\rho^v>^v = \frac{353.4 P_0}{T} - 0.611 <\rho^v_w>^v
 * \f]
 * @param d Oswin isotherm paramters
 * @param X Moisture content [kg/kg db]
 * @param T Temperature [K]
 * @param P0 Pressure [Pa]
 * @returns Density [kg/m^3]
 */
double rho_v(oswin *d, double X, double T, double P0)
{
    double result;

    P0 = P0/101325; /* Convert to atm */

    /* Eq. 3.40 */
    result = 353.4*P0/T - 0.611*rho_v_w(d, X, T);
    return result;
}

double DrhovwDx(oswin *d, double X, double T, double P)
{
    double Pv = vapor_pressure(T), /* Vapor pressure of water */
           aw = OswinInverse(d, X, T), /* Water activity */
           DawDx = OswinDawDx(d, X, T), /* d/dx[a_w] */
           Rw = 461.52, /* (J/kg-K Gas Constant)/(Molar mass of water) */
           result;
    P = P/101325;
    //Pv = Pv/101325;

    /* d/dx[ <rhow>/<rhov> ], solved using Maxima */
    result = (353400000. * Pv * DawDx * P * Rw)
       /(124891560000.*P*P*Rw*Rw - 431854800.*Pv*aw*P*Rw + 373321.*Pv*Pv*aw*aw);
    return result;
}

double permeability(oswin *d, double X, double T)
{
    double eta_w = viscosity(T), /* water viscosity */
           rho_w, /* water density */
           R = GASCONST, /* Gas constant (J/mol-K) */
           Rw = 461.52, /* (J/kg-K Gas Constant)/(Molar mass of water) */
           D0 = 1.78e-5, /* m^2/s */
           //D0 = 149.7e-12,
           Ea = 8.81 * 4184., /* kcal/mol to J/mol */
           Dlfree = D0*exp(-Ea/(R*T)), /* Diffusivity of free water */
           Eb = BindingEnergyOswin(d, X, T), /* Binding energy */
           //n = 2.1; /* Empirical parameter for normal pasta */
           n = 10;
    choi_okos *co;

    co = CreateChoiOkos(WATERCOMP);
    rho_w = rho(co, T);
    DestroyChoiOkos(co);

    /* Eq. 3.53 */
    //return SelfDiffWater(T)*exp(-Eb/(n*R*T)) * eta_w/(rho_w*R*T);
    return Dlfree*exp(-Eb/(n*R*T));// * eta_w/(rho_w*R*T);
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
    double P = 101325, /* 1 atm converted to Pa*/
           K, /* permeability */
           rho_w, /* water density */
           rho_s, /* solid density */
           rhoV, /* <rho^v>^v */
           rhoVW, /* <rho^v_w>^v */
           DrhovwDX, /* d/dx [ <rho^v_w>^v/<rho^v>^v ] */
           epsilon = 0.07, /* Porosity(?) */
           DlnawDX, /* d/dx [ln a_w] */
           eta_w = viscosity(T), /* Water viscosity */
           R = GASCONST, /* Gas constant kg/mol-K */
           Rw = 461.52, /* (J/kg-K Gas Constant)/(Molar mass of water) */
           //Dvap = VaporDiff(T, P), /* Vapor diffusivity from BSL */
           Dvap = 1e-7,
           //Dvap = 5e-8,
           Deff, term1, term2;
    oswin *d;

    d = CreateOswinXiong();

    choi_okos *co;
    co = CreateChoiOkos(PASTACOMP);
    rho_s = rho(co, T);
    DestroyChoiOkos(co);
    co = CreateChoiOkos(WATERCOMP);
    rho_w = rho(co, T);
    DestroyChoiOkos(co);
    //rho_w = 1000;
    //rho_s = 1.47*rho_w;

    //K = 1.5e-16 * exp(-BindingEnergyOswin(d, X, T)/(8.3145*T)),
    K = permeability(d, X, T);
    //K = perm_wat(X*rho_s, epsilon, T);
    rhoV = rho_v(d,X,T,P),
    rhoVW = rho_v_w(d,X,T),
    DrhovwDX = DrhovwDx(d, X, T, P);
    DlnawDX = OswinDlnawDx(d, X, T);


    /* Eq. 3.39 */
    term1 = K*rho_w/(eta_w*(1-epsilon)*rho_s) * (rho_w*Rw*T) * DlnawDX;
    term2 = rhoV*Dvap*DrhovwDX/(rho_w*(1+rhoVW/rhoV));

    printf("P=%g, T=%g, K=%g,\nrhow=%g, rhos = %g,\nrhoV=%g, rhoVW=%g, DrhovwDX=%g,\n"
           "epsilon=%g, DawDx=%g, etaw=%g,\nR=%g, Dvap=%g, Eb=%g\n"
           "term1 = %g, term2 = %g\n\n",
           P,T,K,rho_w,rho_s,rhoV,rhoVW,DrhovwDX,epsilon,
           OswinDawDx(d,X,T),eta_w,R,Dvap,BindingEnergyOswin(d,X,T),term1, term2);
    Deff = term1 + term2;

    DestroyOswinData(d);

    return Deff;
}

double DiffAchantat2(double X, double T)
{
    double P = 101325, /* 1 atm converted to Pa*/
           K, /* permeability */
           rho_w, /* water density */
           rho_s, /* solid density */
           rhoV, /* <rho^v>^v */
           rhoVW, /* <rho^v_w>^v */
           DrhovwDX, /* d/dx [ <rho^v_w>^v/<rho^v>^v ] */
           epsilon = 0.07, /* Porosity(?) */
           DlnawDX, /* d/dx [ln a_w] */
           eta_w = viscosity(T), /* Water viscosity */
           R = GASCONST, /* Gas constant kg/mol-K */
           Rw = 461.52, /* (J/kg-K Gas Constant)/(Molar mass of water) */
           //Dvap = VaporDiff(T, P), /* Vapor diffusivity from BSL */
           Dvap = 1e-7,
           //Dvap = 5e-8,
           Deff, term1, term2;
    oswin *d;

    d = CreateOswinXiong();

    choi_okos *co;
    co = CreateChoiOkos(PASTACOMP);
    rho_s = rho(co, T);
    DestroyChoiOkos(co);
    co = CreateChoiOkos(WATERCOMP);
    rho_w = rho(co, T);
    DestroyChoiOkos(co);
    //rho_w = 1000;
    //rho_s = 1.47*rho_w;

    //K = 1.78e-5 * exp(-BindingEnergyOswin(d, X, T)/(8.314*T)),
    K = permeability(d, X, T);
    //K = perm_wat(X*rho_s, epsilon, T);
    rhoV = rho_v(d,X,T,P),
    rhoVW = rho_v_w(d,X,T),
    DrhovwDX = DrhovwDx(d, X, T, P);
    DlnawDX = OswinDlnawDx(d, X, T);


    /* Eq. 3.39 */
    term1 = K*rho_w/(eta_w*(1-epsilon)*rho_s) * (rho_w*R*T) * DlnawDX;
    term2 = rhoV*Dvap*DrhovwDX/(rho_w*(1+rhoVW/rhoV));

    printf("P=%g, K=%g,\nrhow=%g, rhos = %g,\nrhoV=%g, rhoVW=%g, DrhovwDX=%g,\n"
           "epsilon=%g, DawDx=%g, etaw=%g,\nR=%g, Dvap=%g,\n"
           "term1 = %g, term2 = %g\n\n",
           P,K,rho_w,rho_s,rhoV,rhoVW,DrhovwDX,epsilon,
           OswinDawDx(d,X,T),eta_w,R,Dvap,term1, term2);
    Deff = term1 + term2;

    DestroyOswinData(d);

    return term2;
}
