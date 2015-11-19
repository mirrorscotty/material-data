/**
 * @file binding.c
 * Calculate binding energy based on either the GAB isotherm or the Oswin
 * isotherm. Used in diffusivity calculations.
 */

#include <math.h>

#include "isotherms.h"
#include "constants.h"
#include "diffusivity.h"

#define DX 1e-7

/**
 * Calculate the binding energy based on the GAB equation.
 * @param dat Set of temperature dependent GAB parameters
 * @param X Dry-basis moisture content [kg/kg db]
 * @param T Absolute temperature [K]
 * @returns Binding energy [J/mol]
 *
 * @see BindingEnergyOswin
 */
double BindingEnergyGAB(gab *dat, double X, double T)
{
    double Eb, h, R, dlnawdT;
    h = DX; /* Set the value of dx used for numerical differentiation */
    R = GASCONST; /* Gas constant */

    /* Derivative of ln(aw) with respect to temperature */
    dlnawdT = (log(GABInverse(dat, X, T+h)) - log(GABInverse(dat, X, T-h)))/(2*h);
    /* Multiply by the gas constant and finish differentiating with the chain
     * rule */
    Eb = T*T*R*dlnawdT;

    return Eb;
}

/**
 * Calculate binding energy based on the Oswin isotherm model.
 * Source: Eq 10.21, Handbook of Food Engineering, Second Edition, Ch 10
 * @param dat Set of constants for the isotherm equation.
 * @param X Dry basis moisture content [kg/kg db]
 * @param T Absolute temperature [K]
 * @returns Binding energy [J/mol]
 *
 * @see BindingEnergyGAB
 */
double BindingEnergyOswin(oswin *dat, double X, double T)
{
    double Eb, h, R, dlnawdT;
    h = DX; /* Set the value of dx used for numerical differentiation */
    R = GASCONST; /* Gas constant */

    /* Derivative of ln(aw) with respect to temperature */
    dlnawdT = (log(OswinInverse(dat, X, T+h)) - log(OswinInverse(dat, X, T-h)))/(2*h);

    /* Multiply by the gas constant and finish differentiating with the chain
     * rule */
    Eb = T*T*R*dlnawdT;

    return Eb;
}

/**
 * Calculate binding energy based on the modified Henderson model from
 * Litchfield 1992.
 * @param dat Set of constants for the isotherm equation.
 * @param X Dry basis moisture content [kg/kg db]
 * @param T Absolute temperature [K]
 * @returns Binding energy [J/mol]
 *
 * @see BindingEnergyGAB
 */
double BindingEnergyHenderson(henderson *d, double X, double T)
{
    double Eb, h, R, dlnawdT;
    h = DX; /* Set the value of dx used for numerical differentiation */
    R = GASCONST; /* Gas constant */

    /* Derivative of ln(aw) with respect to temperature */
    dlnawdT = (log(HendersonIsotherm(d, X, T+h)) - log(HendersonIsotherm(d, X, T-h)))/(2*h);

    /* Multiply by the gas constant and finish differentiating with the chain
     * rule */
    Eb = T*T*R*dlnawdT;

    return Eb;
}

/**
 * Calculate the required binding energy to fit the Xiong, et al. diffusivity
 * model to an arbitrary diffusivity.
 */
double BindingEnergyDiff(DiffXiongData *d, double Deff, double T)
{
   double K = d->K,
          Eb,
          R = 8.314; /* Gas Constant */
   double Dself = d->D0*exp(-d->Ea/(R*T)),
          Dr = (Deff/Dself > .99)?.99:(Deff/Dself);

   //Eb = -R*T * log(Dr * K/(1-Dr));
   Eb = -R*T*log(-1*Dr/(K*(Dr-1)));
   return Eb;
}

double BindingEnergyLitchfield(double X, double T)
{
    DiffXiongData *d;
    double Deff = DiffLitchfield(X, T);

    d = CreateDefaultXiongData();
    return BindingEnergyDiff(d, Deff, T);
}

