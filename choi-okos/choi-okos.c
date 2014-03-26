#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "choi-okos.h"
#include "datafile.h"

#ifdef __linux__ /* Not using windows */
#include <errno.h> /* Include the standard error reporting library */
#endif

#define MALLOC_CHECK_ 1

/* TODO:
 * Add functions for convective heat transfer coefficient
 * Plug memory leaks
 * Double check functions to make sure they're accurate
 */

/* Define a variable to hold any error messages returned by the function */
static const char *error = NULL;

choi_okos* CreateChoiOkos(double Mpro, double Mfat, double Mcar,
                          double Mfib, double Mash, double Mwat,
                          double Mice)
{
    choi_okos *co;
    co = (choi_okos*) calloc(sizeof(choi_okos), 1);
    co->Mpro = Mpro;
    co->Mfat = Mfat;
    co->Mcar = Mcar;
    co->Mfib = Mfib;
    co->Mash = Mash;
    co->Mwat = Mwat;
    co->Mice = Mice;

    /* We don't do freezing... yet. */
    /* TODO: Find appropriate values for each of these */
    co->MW_pro = 2.6148e5;
    co->MW_fat = 2e5;
    co->MW_car = 2e5;
    co->MW_fib = 2e5;
    co->MW_ash = 2e5;
    co->MW_wat = 18.01528;

    co->Hfus = 6010;
    co->Tf = 273.15;
    co->R = 8.314;

    return co;
}

void DestroyChoiOkos(choi_okos *co)
{
    free(co);
}

int report_error(const char *str)
{
    #ifdef __linux__
    fprintf(stderr, "%s: %s\n", str, strerror(errno));
    return errno;
    #else
    fprintf(stderr, "Fatal error: %s\n", str);
    return -1;
    #endif
}

/* Functions to actually calculate properties. */

/* ----------------------- Heat Transfer Stuff ----------------------- */

/**
 * Calculate heat capacity using the magic of the Choi-Okos Equations.
 * @param co Composition data
 * @param T Temperature [K]
 * @returns Heat capacity [J/K]
 */
double Cp(choi_okos *co, double T)
{
    double Cp_wat, Cp_ice, Cp_pro, Cp_fat, Cp_car, Cp_ash, Cp_fib;

    T = T-273.15; /* Convert from Kelvin to Celcius */

    if(T >= 0) {
        Cp_wat = 4.1289 + 9.0864e-5*T - 5.4731e-6*pow(T, 2);
    } else {
        Cp_wat = 4.1289 + 5.3062e-3*T - 9.9516e-4*pow(T, 2);
    }
    Cp_ice = 2.0623 + 6.0769e-3*T;
    Cp_pro = 2.0082 + 1.2089e-3*T - 1.3129e-6*pow(T, 2);
    Cp_fat = 1.9842 + 1.4733e-4*T - 4.8008e-6*pow(T, 2);
    Cp_car = 1.5488 + 1.9625e-3*T - 5.9399e-6*pow(T, 2);
    Cp_ash = 1.0926 + 1.8896e-3*T - 3.6817e-6*pow(T, 2);
    Cp_fib = 1.8459 + 1.8306e-3*T - 4.6509e-6*pow(T, 2);
 
    return co->Mpro*Cp_pro + co->Mfat*Cp_fat + co->Mcar*Cp_car
        + co->Mfib*Cp_fib + co->Mash*Cp_ash + co->Mwat*Cp_wat + co->Mice*Cp_ice;
}

/**
 * Calculate the thermal conductivity using the Choi-Okos equations.
 * @param co Composition data
 * @param T Temperature [K] Must be between -40C to 150C
 * @returns Thermal conductivity [W/(m K)]
 */
double k(choi_okos *co, double T)
{
    /* Define all of the local variables needed */
    double k_pro, k_fat, k_car, k_fib, k_ash, k_wat, k_ice;
    double p_pro, p_fat, p_car, p_fib, p_ash, p_wat, p_ice;
    double Xv_pro, Xv_fat, Xv_car, Xv_fib, Xv_ash, Xv_wat, Xv_ice;

    /* Freezing currently borked */
    //double Mi = M_ice(X_ice(T), X_solids());
    double Mi = 0;

    T -= 273.15; /* Convert from Kelvin to Celcius */

    /* Calculate the thermal conductivities of the materials */
    k_pro = 1.7881e-1 + 1.1958e-3*T - 2.7178e-6*pow(T, 2);
    k_fat = 1.8071e-1 - 2.7604e-4*T - 1.7749e-7*pow(T, 2);
    k_car = 2.0141e-1 + 1.3874e-3*T - 4.3312e-6*pow(T, 2);
    k_fib = 1.8331e-1 + 1.2497e-3*T - 3.1683e-6*pow(T, 2);
    k_ash = 3.2962e-1 + 1.4011e-3*T - 2.9069e-6*pow(T, 2);
    k_wat = 5.7109e-1 + 1.762e-3*T - 6.703e-6*pow(T, 2);
    k_ice = 2.2196 - 6.248*10e-3*T + 1.0154e-4*pow(T, 2);

    /* Calculate the densities */
    p_pro = 1.3299e3 - 5.1840e-1*T;
    p_fat = 9.2559e2 - 4.1757e-1*T;
    p_car = 1.5991e3 - 3.1046e-1*T;
    p_fib = 1.3115e3 - 3.6589e-1*T;
    p_ash = 2.4238e3 - 2.8063e-1*T;
    p_wat = 997.18 + 3.1439e-3*T - 3.7574e-3*pow(T, 2);
    p_ice = 916.89 - 1.3071e-1*T;

    /* Determine the volume fraction of each component */
    Xv_pro = (co->Mpro/p_pro) / 
        ((co->Mwat-Mi)/p_wat + Mi/p_ice + co->Mpro/p_pro + co->Mfat/p_fat
         + co->Mcar/p_car + co->Mfib/p_fib + co->Mash/p_ash);
    Xv_fat = (co->Mfat/p_fat) / 
        ((co->Mwat-Mi)/p_wat + Mi/p_ice + co->Mpro/p_pro + co->Mfat/p_fat
         + co->Mcar/p_car + co->Mfib/p_fib + co->Mash/p_ash);
    Xv_car = (co->Mcar/p_car) /
        ((co->Mwat-Mi)/p_wat + Mi/p_ice + co->Mpro/p_pro + co->Mfat/p_fat
         + co->Mcar/p_car + co->Mfib/p_fib + co->Mash/p_ash);
    Xv_fib = (co->Mfib/p_fib) /
        ((co->Mwat-Mi)/p_wat + Mi/p_ice + co->Mpro/p_pro + co->Mfat/p_fat
         + co->Mcar/p_car + co->Mfib/p_fib + co->Mash/p_ash);
    Xv_ash = (co->Mash/p_ash) /
        ((co->Mwat-Mi)/p_wat + Mi/p_ice + co->Mpro/p_pro + co->Mfat/p_fat
         + co->Mcar/p_car + co->Mfib/p_fib + co->Mash/p_ash);
    Xv_wat = ((co->Mwat-Mi)/p_wat) /
        ((co->Mwat-Mi)/p_wat + Mi/p_ice + co->Mpro/p_pro + co->Mfat/p_fat
         + co->Mcar/p_car + co->Mfib/p_fib + co->Mash/p_ash);
    //Xv_ice = (Mi/p_ice) / ((co->Mwat-Mi)/p_wat + Mi/p_ice + co->Mpro/p_pro + co->Mfat/p_fat + co->Mcar/p_car + co->Mfib/p_fib + co->Mash/p_ash);

    /* Calculate the thermal conductivity and return it */
    return k_pro*Xv_pro + k_fat*Xv_fat + k_car*Xv_car + k_fib*Xv_fib
        + k_ash*Xv_ash + k_wat*Xv_wat; //+ k_ice*Xv_ice;
}

/**
 * Calculate density using the Choi-Okos equations.
 * @param co Composition data
 * @param T Temperature [K] Must be between -40C to 150C
 * @param rho has units of [kg/m^3]
 */
double rho(choi_okos *co, double T)
{
    double p_pro, p_fat, p_car, p_fib, p_ash, p_wat, p_ice;

    //double Mi = M_ice(X_ice(T), X_solids());
    double Mi = 0;

    T -= 273.15; /* Convert from Kelvin to Celcius */

    /* Calculate the densities */
    p_pro = 1.3299e3 - 5.1840e-1*T;
    p_fat = 9.2559e2 - 4.1757e-1*T;
    p_car = 1.5991e3 - 3.1046e-1*T;
    p_fib = 1.3115e3 - 3.6589e-1*T;
    p_ash = 2.4238e3 - 2.8063e-1*T;
    p_wat = 997.18 + 3.1439e-3*T - 3.7574e-3*pow(T, 2);
    p_ice = 916.89 - 1.3071e-1*T;

    return 1/(co->Mpro/p_pro + co->Mfat/p_fat + co->Mcar/p_car + co->Mfib/p_fib
            + co->Mash/p_ash + (co->Mwat-Mi)/p_wat + Mi/p_ice);
}

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
