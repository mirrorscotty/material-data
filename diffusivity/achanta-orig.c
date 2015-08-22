#include <math.h>
#include <stdio.h>

double AchantaOrigDeff(double X, double Tk)
{
    double Tc = Tk - 273.15, /* Temperature [C] */
           P = 0.77, /* Pressure [atm] */
           Dv = 1e-10*pow(Tk/328.15, 1.85) * 2/P, /* Vapor diffusivity */
           Pv0 = log(7.967 - 1668.0/(Tk-45.15))/760, /* Vapor pressure */
           eps = 0.27, /* Porosity [-] */
           rhos = 1.21, /* solid density */
           rhorat = 1.21, /* rhos/rhow */
           Rw = 4.56e-3, /* Gas constant [m^3-atm/K-kgH20] */
           R = 1.987, /* Gas constant [kcal/mol-K] */
           n = 10,

           Eb = 24e3*exp(-25*X), /* Binding energy */
           D0 = 149.7e-12,
           Dl = D0 * exp(-Eb/(n*R*Tk)), /* Liquid diffusivity */

           A = Pv0/(Rw*Tk),
           K = 0.176-1.748e-3*Tc, /* Oswin isotherm parameter "K" */
           N = 0.182 + 6.946e-3*Tc, /* Oswin isotherm parameter "N" */
           par1 = pow(K/X, 1/N) + 1,
           par2 = A*(par1-1.0)/(pow(par1,2)*N*X),
           
           rhoVW = A/par1,
           rhoV = 353.4/Tk - 0.611*rhoVW,

           par4 = (353.4/Tk) * par2/pow(rhoV,2),
           H = rhoV*Dv/(1.0e3*(1-eps)*(rhoVW/rhoV)),
           par5 = H*par4/rhorat,
           Deff = Dl + par5;

    return Deff;
}

