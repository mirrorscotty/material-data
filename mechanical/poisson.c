#include <math.h>

/**
 * Calculate the Poisson ratio for an elastic, porous material as a function
 * of porosity and pore-free poisson ratio.
 *
 * This is equation 2 in Arnold, Boccaccini, and Ondracek, 1996, and is
 * appropriate for spherical pore geometries and low porosity.
 *
 * @param phi Porosity [-]
 * @param v0 Poisson ratio at zero porosity [-]
 * @returns Poisson ratio
 */
double poisson(double phi, double v0)
{
    double v;
    v = 0.5-(pow(1-pow(phi, 2./3.), 1.21) * (2*(3-5*phi)*(1-2*v0)+3*phi*(1+v0)))
        / (4*(3-5*phi)*(1-phi));
    return v;
}

