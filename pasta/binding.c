/**
 * @file binding.c
 * Calculate binding energy based on either the GAB isotherm or the Oswin
 * isotherm. Used in diffusivity calculations.
 */

#include <math.h>

#include "isotherms.h"

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
    h = .000001; /* Set the value of dx used for numerical differentiation */
    R = 8.314; /* Gas constant */

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
    h = .000001; /* Set the value of dx used for numerical differentiation */
    R = 8.314; /* Gas constant */

    /* Derivative of ln(aw) with respect to temperature */
    dlnawdT = (log(OswinInverse(dat, X, T+h)) - log(OswinInverse(dat, X, T-h)))/(2*h);

    /* Multiply by the gas constant and finish differentiating with the chain
     * rule */
    Eb = T*T*R*dlnawdT;

    return Eb;
}

