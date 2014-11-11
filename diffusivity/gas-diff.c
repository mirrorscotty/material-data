#include "choi-okos.h"
#include "isotherms.h"
#include "constants.h"
#include "conversions.h"
#include "pasta.h"
#include "diffusivity.h"
#include <math.h>

/**
 * Calculate the diffusivity of water in a nonpolar gas.
 * Source: Transport Phenomena, Second Edition, Eq 17.2-1
 *
 * Calculates the diffusivity at low pressures, and results in approximately a
 * 6-8% error at atmospheric pressure.
 *
 * @param T Temperature [K]
 * @param P Pressure [Pa]
 * @returns Diffusivity [m^2/s]
 */
double VaporDiff(double T, double P)
{
    double D, /* Binary diffusivity [cm^2/s] */
           Tca = 647.14, /* Critical temperature (Water) [K] */
           Tcb = 126.21, /* Critical temperature (Nitrogen) [K] */
           Pca = 217.8, /* Critical pressure (Water) [atm] */
           Pcb = 33.46, /* Critical pressure (Nitrogen) [atm] */
           Ma = 18.0153, /* Molar mass (Water) [g/mol] */
           Mb = 28.0134, /* Molar mass (Nitrogen) [g/mol] */
           a = 3.64e-4, /* Dimensionless parameter for H2O with a nonpolar gas*/
           b = 2.334; /* Dimensionless parameter for H2O with a nonpolar gas */

    P = PA_TO_ATM(P); /* Convert from Pa to atm */

    /* The temperatures for this equation must be in Kelvin, and the pressures
     * should be in Pascals. The diffusivity is calculated in cm^2/s */
    D = a*pow(T/sqrt(Tca*Tcb), b);
    D *= (pow(Pca*Pcb, 1/3.) * pow(Tca*Tcb, 5/12.) * sqrt(1/Ma+1/Mb))/P;

    D *= 1e-4; /* Convert to m^2/s */

    return D;
}


