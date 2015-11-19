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
    //result = 1./760. * log(7.967 - 1668.0/(T-45.15));
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

//double permeability(oswin *d, double X, double T)
//{
//    double eta_w = viscosity(T), /* water viscosity */
//           rho_w, /* water density */
//           R = GASCONST, /* Gas constant (J/mol-K) */
//           Rw = 461.52, /* (J/kg-K Gas Constant)/(Molar mass of water) */
//           //D0 = 149.7e-12,
//           D0 = 1.78e-5, /* m^2/s */
//           Ea = 8.81 * 4184., /* kcal/mol to J/mol */
//           Dlfree = D0*exp(-Ea/(R*T)), /* Diffusivity of free water */
//           Eb = BindingEnergyOswin(d, X, T), /* Binding energy */
//           //n = 2.1; /* Empirical parameter for normal pasta */
//           n = 10;
//    choi_okos *co;
//
//    co = CreateChoiOkos(WATERCOMP);
//    rho_w = rho(co, T);
//    DestroyChoiOkos(co);
//
//    /* Eq. 3.53 */
//    //return SelfDiffWater(T)*exp(-Eb/(n*R*T)) * eta_w/(rho_w*R*T);
//    return Dlfree*exp(-Eb/(n*R*T));// * eta_w/(rho_w*R*T);
//}

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
           rho_w, /* water density */
           rho_s, /* solid density */
           rhoV, /* <rho^v>^v */
           rhoVW, /* <rho^v_w>^v */
           DrhovwDX, /* d/dx [ <rho^v_w>^v/<rho^v>^v ] */
           epsilon = 0.07, /* Porosity(?) */
           R = GASCONST, /* Gas constant kg/mol-K */
           Rw = 461.52, /* (J/kg-K Gas Constant)/(Molar mass of water) */
           D0 = 1.78e-5, /* m^2/s */
           Ea = 8.81 * 4184., /* kcal/mol to J/mol */
           Eb,
           n = 2.1,
           Dlfree = D0*exp(-Ea/(R*T)), /* Diffusivity of free water */
           //Dvap = VaporDiff(T, P), /* Vapor diffusivity from BSL */
           //Dvap = 1e-7,
           Dvap = 1e-10*pow(T/238.15, 1.85) * 2/(P/101325),
           //Dvap = 5e-8,
           Deff, term1, term2;
    oswin *d;

    d = CreateOswinXiong();

    /* Calculate the densities for both water and dry pasta at the desired
     * temperature */
    choi_okos *co;
    co = CreateChoiOkos(PASTACOMP);
    rho_s = rho(co, T);
    DestroyChoiOkos(co);
    co = CreateChoiOkos(WATERCOMP);
    rho_w = rho(co, T);
    DestroyChoiOkos(co);
    //rho_w = 1000;
    //rho_s = 1.47*rho_w;

    Eb = BindingEnergyOswin(d, X, T), /* Binding energy */
    rhoV = rho_v(d,X,T,P),
    rhoVW = rho_v_w(d,X,T),
    DrhovwDX = DrhovwDx(d, X, T, P);

    /* Eq. 3.39 */
    //term1 = K*rho_w/(eta_w*(1-epsilon)*rho_s) * (rho_w*Rw*T) * DlnawDX;
    term1 = Dlfree*exp(-Eb/(n*R*T));
    term2 = rhoV*Dvap*DrhovwDX/(rho_w*(1-epsilon)*(rhoVW/rhoV));

    Deff = term1 + term2;

    DestroyOswinData(d);

    return Deff;
}

/**
 * Function to give to fitnlmM to calculate the D0, Ea, and Dvap values that
 * best fit a set of data using nonlinear regression.
 * @param x Row matrix of independent values for the function. x[0,0]: Moisture
 *      Content [kg/kg db] x[0,1]: Temperature [K]
 * @param beta Column matrix of the square root of fitting parameters. B[0,0]:
 *      D0 B[1,0]: Ea B[2,0]:Dvap
 * @returns Effective diffusivity [m^2/s]
 */
double AchantaDiffModel(matrix *x, matrix *beta)
{
    double D0 = val(beta, 0, 0), /* m^2/s */
           Ea = val(beta, 1, 0), /* J/mol */
           Dvap = val(beta, 2, 0); /* m^2/s */
    double X = val(x, 0, 0), /* kg/kg db */
           T = val(x, 0, 1) + 273.15; /* K */

    /* Square the values we're given. The program calling this function passes
     * the square root of these values, and this ensures that they remain
     * positive */
    D0 = pow(D0,2);
    Ea = pow(Ea,2);
    Dvap = pow(Dvap,2);

    double P = 101325, /* 1 atm converted to Pa*/
           rho_w, /* water density */
           rho_s, /* solid density */
           rhoV, /* <rho^v>^v */
           rhoVW, /* <rho^v_w>^v */
           DrhovwDX, /* d/dx [ <rho^v_w>^v/<rho^v>^v ] */
           epsilon = 0.07, /* Porosity(?) */
           R = GASCONST, /* Gas constant kg/mol-K */
           Rw = 461.52, /* (J/kg-K Gas Constant)/(Molar mass of water) */
           Eb,
           n = 2.1,
           Dlfree = D0*exp(-Ea/(R*T)), /* Diffusivity of free water */
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

    Eb = BindingEnergyOswin(d, X, T), /* Binding energy */
    rhoV = rho_v(d,X,T,P),
    rhoVW = rho_v_w(d,X,T),
    DrhovwDX = DrhovwDx(d, X, T, P);

    term1 = Dlfree*exp(-Eb/(n*R*T));
    term2 = rhoV*Dvap*DrhovwDX/(rho_w*(1-epsilon)*(rhoVW/rhoV));

    Deff = term1 + term2;

    DestroyOswinData(d);

    return Deff;
}

double DiffAchantaFitted(double X, double T)
{
    double Deff;
    matrix *x, *beta;

    x = CreateMatrix(1, 2);
    setval(x, X, 0, 0);
    setval(x, T, 0, 1);

    beta = CreateMatrix(3, 1);
    setval(beta, 0.000171135, 0, 0);
    setval(beta, 41640, 1, 0);
    setval(beta, 1.32835e-09, 2, 0);

    Deff = AchantaDiffModel(x, beta);

    DestroyMatrix(x);
    DestroyMatrix(beta);

    return Deff;
}

