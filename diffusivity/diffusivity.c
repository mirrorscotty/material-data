/**
 * @file diffusivity.c
 * Several models for calculating diffusivity in porous media.
 */

#include "choi-okos.h"
#include "isotherms.h"
#include "constants.h"
#include "conversions.h"
#include "pasta.h"
#include "diffusivity.h"
#include <math.h>

/** 
 * Calculate diffusivity in pasta based on the model outlined in chapter 10 of
 * the Handbook of Food Engineering, Second Edition. This function uses the
 * Oswin isotherm to determine binding energy.
 * @param X Moisture Content [kg/kg db]
 * @param T Temperature [K]
 * @returns Diffusivity [m^2/s]
 *
 * @see DiffCh10GAB
 */
double DiffCh10(double X, double T)
{
    oswin *dat;
    dat = OSWINDATA();

    double Deff,
           D0 = 6.3910e-8, /* [m^2/s] Source: Handbook of Food Engineering */
           //Ea = 25900, /* Source: Litchfield and Okos (1986) */
           Ea = 21760, /* [J/mol] Source: Xiong et al. (1991) */
           K = 1032.558, /* Source: Xiong et al. (1991) */
           Eb = BindingEnergyOswin(dat, X, T),
           R = 8.314; /* Gas Constant */

    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    Deff = D0 * exp(-Ea/(R*T))
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

double DiffCh10new(oswin *dat, double X, double T)
{
    double Deff,
           D0 = 6.3910e-8, /* [m^2/s] Source: Handbook of Food Engineering */
           //Ea = 25900, /* Source: Litchfield and Okos (1986) */
           Ea = 21760, /* [J/mol] Source: Xiong et al. (1991) */
           K = 1032.558, /* Source: Xiong et al. (1991) */
           Eb = BindingEnergyOswin(dat, X, T),
           R = 8.314; /* Gas Constant */

    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    Deff = D0 * exp(-Ea/(R*T))
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

/**
 * Modification of the diffusivity equation from Xiong et al. 1991 to use
 * the self-diffusivity of water.
 */
double DiffCh10Mod(double X, double T)
{
    oswin *dat;
    dat = OSWINDATA();

    double Deff,
           Dself = SelfDiffWater(T),
           phi = POROSITY,
           tau = 10.5409, /* Fitted value for tortuosity based on 55 C data */
           Sw = sat_wat(conc_wat(X, POROSITY, T), POROSITY, T),
           K = 1032.558, /* Source: Xiong et al. (1991) */
           Eb = BindingEnergyOswin(dat, X, T),
           R = 8.314; /* Gas Constant */

    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    Deff = phi/tau * Dself
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

/** 
 * Calculate diffusivity in pasta based on the model outlined in chapter 10 of
 * the Handbook of Food Engineering, Second Edition. This function uses the
 * GAB isotherm to determine binding energy.
 * @param X Moisture Content [kg/kg db]
 * @param T Temperature [K]
 * @returns Diffusivity [m^2/s]
 *
 * @see DiffCh10
 */
double DiffCh10GAB(double X, double T)
{
    gab *dat;
    dat = GABDATA();

    double Deff,
           D0 = 6.3910e-8, /* Source: Xiong et al (1991) */
           //Ea = 25900, /* Source: Litchfield and Okos (1986) */
           Ea = 21760, /* [J/mol] Source: Xiong et al. (1991) */
           K = 1032.6, /* Source: Xiong et al. (1991) */
           Eb = BindingEnergyGAB(dat, X, T),
           R = 8.314; /* Gas Constant */


    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    Deff = D0 * exp(-Ea/(R*T))
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

/* Calculate diffusivity based on the model in Chapter 10 of The Handbook of
 * Food Engineering. The isotherm model used is the modified Henderson isotherm
 * described in Litchfield 1992.
 * @param X Moisture Content [kg/kg db]
 * @param T Temperature [K]
 * @returns Diffusivity [m^2/s]
 *
 * @see DiffCh10
 */
double DiffCh10Hend(double X, double T)
{
    henderson *dat;
    dat = CreateHendersonData();

    double Deff,
           D0 = 6.3910e-8, /* Source: Xiong et al (1991) */
           //Ea = 25900, /* Source: Litchfield and Okos (1986) */
           Ea = 21760, /* [J/mol] Source: Xiong et al. (1991) */
           K = 1032.6, /* Source: Xiong et al. (1991) */
           Eb = BindingEnergyHenderson(dat, X, T),
           R = GASCONST; /* Gas Constant */


    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    Deff = D0 * exp(-Ea/(R*T))
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

/**
 * Calculate vapor diffusivity based on the Oswin isotherm and the Chapter 10
 * model. The only difference is that the D0 value is somewhat larger for vapor
 * diffusivity than it is for liquid diffusivity.
 * @param X Moisture content [kg/kg db]
 * @param T Temperature [K]
 * @returns Diffusivity [m^2/s]
 *
 * @see DiffCh10
 */
double VaporDiffCh10(double X, double T)
{
    double D0gas = 2e-5,
           D0liq = 6.3910e-8;
    return D0gas/D0liq * DiffCh10(X, T);
}

/**
 * Calculate the diffusivity for unsaturated porous materials using the model
 * found in Xiong et. al.
 * @param X Moisture content [kg/kg db]
 * @param T Temperature [K]
 * @param phi Porosity [-]
 * @returns Effective diffusion coefficient [m^2/s]
 */
double PorousDiffCh10(double X, double T, double phi)
{
    double D0gas = 2e-5,
           D0liq = 6.3910e-8,
           D;
    D = DiffCh10(X, T);
    return (phi * D0gas/D0liq + (1-phi)) * D;
}

double PorousDiff(double X, double T, double phi)
{
    double P = 101300, /* [Pa] */
           Dliq, Dgas;
    Dliq = DiffCh10(X, T);
    Dgas = VaporDiff(T, P);
    return phi*Dgas + (1-phi)*Dliq;
}

/**
 * Temperature-dependent self diffusion coefficent for water.
 * Valid from 0C to 100C. From Holz et al. 2000
 * @param T Temperature [K]
 * @returns Diffusion coefficient [m^2/s]
 */
double SelfDiffWater(double T)
{
    double D,
           D0 = 1.635e-8,
           Ts = 215.05,
           gamma = 2.063;
    D = D0*pow(T/Ts-1, gamma);
    return D;
}

double KnudsenDiff(double T)
{
    double D,
           a = 0, /* Pore radius */
           R = GASCONST, /* Gas constant */
           pi = M_PI, /* Value of pi */
           Mw = 18.0153; /* Molar mass of water (g/mol) */

    D = 8*a/3 * sqrt(R*T/(2*pi*Mw));

    return D;
}

/**
 * Calculate diffusivity from length and the kF value determined from the Crank
 * equation.
 * @param kF (pi^2 D)/L^2 [1/s]
 * @param L Sample thickness [m]
 * @returns Diffusivity [m^2/s]
 */
double DiffkF(double kF, double L)
{
    double D;

    D = kF*L*L/(M_PI*M_PI);
    
    return D;
}

