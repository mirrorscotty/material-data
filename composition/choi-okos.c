/**
 * @file choi-okos.c
 * Implementation of the Choi-Okos equations for calculating food properties
 * based on composition.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "choi-okos.h"

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
    /* Define the proximate composition of the material */
    co->Mpro = Mpro;
    co->Mfat = Mfat;
    co->Mcar = Mcar;
    co->Mfib = Mfib;
    co->Mash = Mash;
    co->Mwat = Mwat;
    co->Mice = Mice;

    /* We don't do freezing... yet. */
    /* TODO: Find better values. Current ones are from Choi-Okos (1986). */
    co->MW_pro = 3e4; /* Whey protein */
    co->MW_fat = 3e4; /* Assume this is the same as protein for no reason */
    co->MW_car = 3e5; /* Starch */
    co->MW_fib = 3e5; /* Cellulose */
    co->MW_ash = 158; /* Milksalt */
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

/**
 * Take a set of composition data for dry solids and recalculate the mass
 * fractions of each component to account for water.
 * @param co Dry composition
 * @param Xdb Dry basis moisture content [kg/kg db]
 * @return Wet composition
 */
choi_okos* AddDryBasis(choi_okos *co, double Xdb)
{
    choi_okos *cowet;
    cowet = CreateChoiOkos(0, 0, 0, 0, 0, 0, 0);
    /* Assume that there's no water or ice initially. */
    cowet->Mpro = co->Mpro/(1+Xdb);
    cowet->Mfat = co->Mfat/(1+Xdb);
    cowet->Mcar = co->Mcar/(1+Xdb);
    cowet->Mfib = co->Mfib/(1+Xdb);
    cowet->Mash = co->Mash/(1+Xdb);
    cowet->Mwat = Xdb/(1+Xdb);

    return cowet;
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

