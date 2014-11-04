/**
 * @file freezing.c
 * Functions for calculating food properties as the moisture in the samples is
 * in the process of freezing.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "choi-okos.h"

/* ----------------------- Freezing Stuff ----------------------- */
/* If I'm not lazy and I've actually fixed these functions, compile them. */
#ifdef NOT_LAZY

/**
 * Calculate the mole fraction of "a" given it's mass fraction and molecular
 * weight. The total number of mole is calculated from the composition data
 * provided in the data file.
 */
double MoleFrac(double Ma, double MWa)
{
    double total_moles;
    total_moles = (Mwat+Mice)/MW_wat + Mpro/MW_pro + Mcar/MW_car + Mfib/MW_fib + Mash/MW_ash;
    return (Ma/MWa)/total_moles;
}

/**
 * Determine the average molecular weight of the soluble solids. Fat is not
 * soluble.
 */
double MW_solids()
{
    double  Xpro = MoleFrac(Mpro, MW_pro),
            Xcar = MoleFrac(Mcar, MW_car),
            Xfib = MoleFrac(Mfib, MW_fib),
            Xash = MoleFrac(Mash, MW_ash),
            Xs = Xpro + Xcar + Xfib + Xash;


    return (MW_pro*Xpro + MW_car*Xcar + MW_fib*Xfib + MW_ash*Xash)/Xs;
}

/**
 * Return the average mass fraction of the aqueous phase.
 */
double MW_average()
{
    double  Xpro = MoleFrac(Mpro, MW_pro),
            Xcar = MoleFrac(Mcar, MW_car),
            Xfib = MoleFrac(Mfib, MW_fib),
            Xash = MoleFrac(Mash, MW_ash),
            Xwat = MoleFrac(Mwat, MW_wat);
    return Xwat*MW_wat + Xpro*MW_pro + Xcar*MW_car + Xfib*MW_fib + Xash*MW_ash;
}

/* Only have these functions work if we're using the freezing model. Adding
 * this stuff to the sterilization model kind of breaks it. */
#ifdef CALC_ICE_FORMATION
/**
 * Calculate the mole fraction of ice in food given the temperature and mole
 * fraction of solids.
 * T has units of Kelvin.
 */
double X_ice(double T)
{
    /* Calculate the initial freezing temperature */
    double Ti, x_w1, Xs;

    Xs = X_solids();
    Ti = pow( (1/Tf - R/Hfus*log(1-Xs)), -1);

    if(T>Ti) {
        return 0; /* No ice formed above the freezing point */
    } else {
        x_w1 = exp((1/Tf - 1/T) * Hfus/R);
        return ( 1-x_w1-Xs );
    }
}

/**
 * @brief IceMassFrac
 * @param T Temperature (in K)
 * @return Mass fraction of ice.
 */
double IceMassFrac(double T) {
    return X_ice(T)*MW_wat/MW_average();
}

/**
 * @brief Calculate the mole fraction of soluble solids.
 * @return Mole fraction of solids.
 */
double X_solids()
{
    double  Xpro = MoleFrac(Mpro, MW_pro),
            Xcar = MoleFrac(Mcar, MW_car),
            Xfib = MoleFrac(Mfib, MW_fib),
            Xash = MoleFrac(Mash, MW_ash),
            Xs = Xpro + Xcar + Xfib + Xash;
    return Xs;
}

/**
 * @brief Calculate the mass fraction of ice given the mole fraction of ice and solids.
 * @param x Mole fraction of ice/water
 * @param y Mole fraction of solids
 * @return Mass fraction of ice.
 */
double M_ice(double x, double y)
{
    double MW_s = MW_solids();
    return ( x*MW_wat/(x*MW_wat+(1-x-y)*MW_wat+y*MW_s) );
}

#else /* Return safe values if we're not doing freezing */
double X_ice(double T) { return 0; }
double X_solids() { return 1; }
double M_ice(double x, double y) { return 0; }
#endif

/* Return the density of all of the solids (excluding ice, which is calculated
 * separately. */
/**
 * @brief Determine the density of all the solids.
 * Ice is excluded from this calculation and is determined separately.
 * @param T Temperature (in K)
 * @return Average density (kg/m^3)
 */
double p_solids(double T)
{
    double p_pro, p_fat, p_car, p_fib, p_ash;
    T = T-273.15;

    /* Calculate the densities */
    p_pro = 1.3299e3 - 5.1840e-1*T;
    p_fat = 9.2559e2 - 4.1757e-1*T;
    p_car = 1.5991e3 - 3.1046e-1*T;
    p_fib = 1.3115e3 - 3.6589e-1*T;
    p_ash = 2.4238e3 - 2.8063e-1*T;

    return 1/(Mpro/p_pro + Mfat/p_fat + Mcar/p_car + Mfib/p_fib + Mash/p_ash);
}

/* Density of liquid water. */
double p_water(double T)
{
    T = T-273.15;
    return 997.18 + 3.1439e-3*T - 3.7574e-3*pow(T, 2);
}

/* Density of ice. */
double p_ice(double T)
{
    T = T-273.15;
    return 916.89 - 1.3071e-1*T;
}

/**
 * @brief Determine the volume fraction of water.
 * @param T Temperature (K)
 * @return Volume fraction of water
 */
double Xv_water(double T)
{
    double Mi, Mw, Ms;

    Mi = M_ice(X_ice(T), X_solids());
    Mw = M_ice(MoleFrac(Mwat, MW_wat)-X_ice(T), X_solids());
    Ms = 1-Mi-Mw;

    /* FixMe! */
    return (Mw/p_water(T)) / (Mi/p_ice(T) + Mw/p_water(T) + Ms/p_solids(T));
}

/**
 * @brief Determine the volume fraction of ice.
 * @param T Temperature (K)
 * @return Volume fraction of ice
 */
double Xv_ice(double T)
{
    double Mi, Mw, Ms;

    Mi = M_ice(X_ice(T), X_solids());
    Mw = M_ice(MoleFrac(Mwat, MW_wat)-X_ice(T), X_solids());
    Ms = 1-Mi-Mw;

    /* FixMe! */
    return (Mi/p_ice(T)) / (Mi/p_ice(T) + Mw/p_water(T) + Ms/p_solids(T));
}

/**
 * @brief Calculate thermal diffusivity during freezing.
 * @param T Temperature (K)
 * @return Thermal diffusivity (m^2/s)
 */
double alphaFZ(double T)
{
    return k(T)/(rho(T)*CpFz(T));
}

/**
 * @brief Calculate the heat capacity of the partially frozen food products.
 * TODO: Compare this against the regular Cp function for temperatures above
 * freezing.
 * @param T Temperature (K)
 * @return Heat capacity ( W/(m K) )
 */
double CpFz(double T)
{
    double Mi, Mw, Ms, dMi, dMw, Ti;
    double Xw, Xs;
    //double dT = 0.0001; /* Used for calculating derivatives */
    double dT = 1;

    Xw = MoleFrac(Mwat, MW_wat);
    Xs = X_solids();

    /* Calculate the initial freezing temperature */
    Ti = pow( (1/Tf - R/Hfus*log(1-Xs)), -1);

    Mi = M_ice(X_ice(T), Xs); /* Mass fraction of ice */
    Mw = M_ice(Xw-X_ice(T), Xs); /* Mass fraction of water */
    Ms = 1-Mi-Mw; /* Mass fraction of solids */

    dMi = (M_ice(X_ice(T+dT), Xs)-Mi)/dT;
    dMw = (M_ice(Xw-X_ice(T+dT), Xs)-Mw)/dT;

    return ( Mw*Cp_water(T) + Ms*Cp_solids(T) + Mi*Cp_ice(T)  - Hfus*dMi -
             (dMw*Cp_water(T) + dMi*Cp_ice(T))*(Ti-T) );
}

/**
 * @brief Heat capacity of just liquid water.
 * @param T Temperature (K)
 * @return Heat capacity ( W/(m K) )
 */
double Cp_water(double T)
{
    T = T-273.15;
    if(T >= 0) {
        return 4.1289 + 9.0864e-5*T - 5.4731e-6*pow(T, 2);
    } else {
        return 4.1289 + 5.3062e-3*T - 9.9516e-4*pow(T, 2);
    }
}

/**
 * @brief Heat capacity of all solids, excluding ice.
 * @param T Temperature (K)
 * @return Heat capacity ( W/(m K) )
 */
double Cp_solids(double T)
{
    double Cp_pro, Cp_fat, Cp_car, Cp_ash, Cp_fib;
    T = T-273.15;
    Cp_pro = 2.0082 + 1.2089e-3*T - 1.3129e-6*pow(T, 2);
    Cp_fat = 1.9842 + 1.4733e-4*T - 4.8008e-6*pow(T, 2);
    Cp_car = 1.5488 + 1.9625e-3*T - 5.9399e-6*pow(T, 2);
    Cp_ash = 1.0926 + 1.8896e-3*T - 3.6817e-6*pow(T, 2);
    Cp_fib = 1.8459 + 1.8306e-3*T - 4.6509e-6*pow(T, 2);

    return Mpro*Cp_pro + Mfat*Cp_fat + Mcar*Cp_car + Mfib*Cp_fib + Mash*Cp_ash;
}

/**
 * @brief Heat capacity of ice
 * @param T Temperature (K)
 * @return Heat capacity ( W/(m K) )
 */
double Cp_ice(double T)
{
    T = T-273.15;
    return 2.0623 + 6.0769e-3*T;
}
#endif

