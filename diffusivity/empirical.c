#include <math.h>

/**
 * Diffusivity model for potato from Chemkhi et al. (2004)
 *
 * Other potato drying papers:
 * Rahman 2003
 * 
 */
double DiffChemkhi(double X, double T)
{
    double a = 1.29e-6,
           b = 0.0725,
           c = 2044;
    return a*exp(-b/X)*exp(-c/T);
}

/**
 * Diffusivity model for pasta from Litchfield and Okos (1992)
 */
double DiffLitchfield(double X, double T)
{
    double A = 2.3920e5,
           B = 3.1563e3, /* Activation energy divided by gas constant */
           C = 7.9082e14,
           D = 1.5706e1,
           E = 6.8589e-1;
    //T=T-273.15; /* Convert from Kelvin to Celcius */
    return A*exp(-B/T)*(1-exp(-C*pow(X, D))+pow(X, E));
}

