#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "can.h"
#include "datafile.h"

#ifdef __linux__ /* Not using windows */
#include<errno.h> /* Include the standard error reporting library */
#endif

#define MALLOC_CHECK_ 1

/* TODO:
 * Add functions for convective heat transfer coefficient
 * Plug memory leaks
 * Double check functions to make sure they're accurate
 */


/* Define a variable to hold any error messages returned by the function */
static const char *error = NULL;

/* Annoying global variable definitions. For now, all the model parameters are
 * hard coded into the dll using these variables. This should be changed to
 * allow them to be read in from a data file.
 */
double Mpro, Mfat, Mcar, Mfib, Mash, Mwat, Mice;
double AA, EaA, AB, EaB, R;
double A, B, C;
double Pressure, molar_mass;
double Sutherland, Tref, muref;
double To, Text_hot, Text_cold;
double v, L, t_heat;

/* Variables for the finite difference solver */
double Deltax, NNodes, Deltat, NTimeSteps;

//int main(int argc, char *argv[])
//{
//    initialize_variables();
//    init("can_data.dat");
	/*output_data()*/
    //printf("%f\n", k(80+273)/(rho(80+273)*Cp(80+273)));
//    printf("%f\n", k(80+273));
//	return 0;
//}

/* Set all global variables to an initial value of zero in case something goes
 * horrible, horribly wrong and something tries to read an uninitialized value.
 */
void initialize_variables()
{
    Mpro = 0;
    Mfat = 0;
    Mcar = 0;
    Mfib = 0;
    Mash = 0;
    Mwat = 0;
    Mice = 0;
    AA = 0;
    EaA = 0;
    AB = 0;
    EaB = 0;
    R = 0;
    A = 0;
    B = 0;
    C = 0;
    Pressure = 0;
    molar_mass = 0;
    Sutherland = 0;
    Tref = 0;
    muref = 0;
    To = 0;
    Text_hot = 0;
    Text_cold = 0;
    v = 0;
    L = 0;
    t_heat = 0;

    Deltax = 0;
    NNodes = 0;
    Deltat = 0;
    NTimeSteps = 0;
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

/* Test function to spit out a table of data with the x values in one column
 * and the results of a function the other.
 */

int output_data()
{
	double min, max;
	int points, i;
	const char *file = "out.csv";

	FILE *fp;
	double x, cp_, k_, rho_;

    min = 0;
    max = 0;
    points = 0;
    i = 0;
    fp = NULL;
    x = 0;
    cp_ = 0;
    k_ = 0;
    rho_ = 0;

	fp = fopen(file, "w");
    if(!fp) {
        report_error("Failed to open file for writing.");
    }

	min = 300;
	max = 480;
	points = 100;

	for (i=1; i <= points; i++) {
		x = min+(max-min)*i/points;
		cp_ = Cp(x);
        k_ = k(x);
        rho_ = rho(x);
		fprintf(fp, "%f,%f,%f,%f\n", x, cp_, k_, rho_);
	}

    if(fp) {
        if(fclose(fp) != 0) {
           report_error("Failed to close file.");
           exit(1);           
        }
    }

	return 0;
}

/* The following three functions are required in order to interface correctly
 * with Comsol
 */
/**
 * Initialization function for the library.
 */
EXTCAN_API int init(const char *str)
{
    return get_vars(str);
}

/**
 * Error reporting function. Returns the last error message set.
 */
EXTCAN_API const char * getLastError()
{
	return error;
}

/**
 * Function to that actually does the work for the library. Basically, Comsol
 * calls this function and passes the name of the function it actually wants
 * to run as an argument. The eval function is responsible for calling the
 * correct function from the library.
 *
 * The LOOP macro simply calls the function multiple times since Comsol expects
 * the result to be an array with the same number of elements 
 */
/*
EXTCAN_API int eval(const char *func,
		   int nArgs,
		   const double **inReal,
		   const double **inImag,
		   int blockSize,
		   double *outReal,
		   double *outImag)
{
	int i;
	if(strcmp(func, "Cp") == 0) {
		LOOP(Cp)
	} else if(strcmp(func, "rho") == 0) {
		LOOP(rho)
	} else if(strcmp(func, "reaction_rate1") == 0) {
		LOOP2(reaction_rate1)
	} else if(strcmp(func, "reaction_rate2") == 0) {
        LOOP2(reaction_rate2)
    } else if(strcmp(func, "k") == 0) {
		LOOP(k)
	} else if(strcmp(func, "T_ext") == 0) {
        LOOP(T_ext)
    } else if(strcmp(func, "h") == 0) {
        LOOP(h)
    } else if(strcmp(func, "T_init") == 0) {
        LOOP(T_init)
    } else if(strcmp(func, "mu") == 0) {
        LOOP(mu)
    } else {
		error = "Cannot find function";
		return 0;
	}
    return 1;
}
*/

/* Functions to actually calculate stuff. */

/**
 * Calculate heat capacity using the magic of the Choi-Okos Equations.
 * T is in Kelvins
 */
double Cp(double T)
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
 
    return Mpro*Cp_pro + Mfat*Cp_fat + Mcar*Cp_car + Mfib*Cp_fib + Mash*Cp_ash + Mwat*Cp_wat + Mice*Cp_ice;
}

/* Calculate the thermal conductivity using the Choi-Okos equations. */
double k(double T)
{
    /* Define all of the local variables needed */
    double k_pro, k_fat, k_car, k_fib, k_ash, k_wat, k_ice;
    double p_pro, p_fat, p_car, p_fib, p_ash, p_wat, p_ice;
    double Xv_pro, Xv_fat, Xv_car, Xv_fib, Xv_ash, Xv_wat, Xv_ice;

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
    Xv_pro = (Mpro/p_pro) / (Mwat/p_wat+Mice/p_ice+Mpro/p_pro+Mfat/p_fat+Mcar/p_car+Mfib/p_fib+Mash/p_ash);
    Xv_fat = (Mfat/p_fat) / (Mwat/p_wat+Mice/p_ice+Mpro/p_pro+Mfat/p_fat+Mcar/p_car+Mfib/p_fib+Mash/p_ash);
    Xv_car = (Mcar/p_car) / (Mwat/p_wat+Mice/p_ice+Mpro/p_pro+Mfat/p_fat+Mcar/p_car+Mfib/p_fib+Mash/p_ash);
    Xv_fib = (Mfib/p_fib) / (Mwat/p_wat+Mice/p_ice+Mpro/p_pro+Mfat/p_fat+Mcar/p_car+Mfib/p_fib+Mash/p_ash);
    Xv_ash = (Mash/p_ash) / (Mwat/p_wat+Mice/p_ice+Mpro/p_pro+Mfat/p_fat+Mcar/p_car+Mfib/p_fib+Mash/p_ash);
    Xv_wat = (Mwat/p_wat) / (Mwat/p_wat+Mice/p_ice+Mpro/p_pro+Mfat/p_fat+Mcar/p_car+Mfib/p_fib+Mash/p_ash);
    Xv_ice = (Mice/p_ice) / (Mwat/p_wat+Mice/p_ice+Mpro/p_pro+Mfat/p_fat+Mcar/p_car+Mfib/p_fib+Mash/p_ash);

    /* Calculate the thermal conductivity and return it */
    return k_pro*Xv_pro + k_fat*Xv_fat + k_car*Xv_car + k_fib*Xv_fib + k_ash*Xv_ash + k_wat*Xv_wat + k_ice*Xv_ice;
}

/**
 * Calculate the reaction rate constant.
 */
double reaction_rate1(double T, double c)
{
    // Fix since I'm too lazy to load it from the GUI;
    R = 8.314;
	return AA*exp(-EaA/(R*T));
}

double reaction_rate2(double T, double c)
{
    // Fix since I'm too lazy to load it from the GUI;
    R = 8.314;
	return AB*exp(-EaB/(R*T));
}

/* Calculate density using the Choi-Okos equations. */
double rho(double T)
{
    double p_pro, p_fat, p_car, p_fib, p_ash, p_wat, p_ice;

    /* Calculate the densities */
    p_pro = 1.3299e3 - 5.1840e-1*T;
    p_fat = 9.2559e2 - 4.1757e-1*T;
    p_car = 1.5991e3 - 3.1046e-1*T;
    p_fib = 1.3115e3 - 3.6589e-1*T;
    p_ash = 2.4238e3 - 2.8063e-1*T;
    p_wat = 997.18 + 3.1439e-3*T - 3.7574e-3*pow(T, 2);
    p_ice = 916.89 - 1.3071e-1*T;

    return 1/(Mpro/p_pro + Mfat/p_fat + Mcar/p_car + Mfib/p_fib + Mash/p_ash + Mwat/p_wat + Mice/p_ice);
}

/* Function to return the external temperature of the can. Based on the time,
 * this is either the hot or cold temperature, so that the can is able to be
 * heated and cooled in a single simulation run.
 */
double T_ext(double t)
{
    if(t <= t_heat) {
        return Text_hot;
    } else {
        return Text_cold;
    }
}

/* Return the initial temperature. Used for setting up the temperature inside
 * the can at the start of the simulation.
 */
double T_init(double t)
{
    return To;
}

/* Calculate the viscosity using an Arrhenius-type equation. */
double mu(double T)
{
    return A*pow(10, (B/(T-C)));
    /* Source: http://en.wikipedia.org/wiki/Viscosity/Viscosity_of_water */
}

/* Determine the convective heat transfer coefficient */
//double h(double T)
//{
//    double Re, Pr, p, Cp_wat, k;
//   T = T-273.15; /* Convert to Celcius */
//
//    /* Calculate the required data for water */
//   k = 5.7109e-1 + 1.762e-3*T - 6.703e-6*pow(T, 2);
//    if(T >= 0) {
//        Cp_wat = 4.1289 + 9.0864e-5*T - 5.4731e-6*pow(T, 2);
//    } else {
//        Cp_wat = 4.1289 + 5.3062e-3*T - 9.9516e-4*pow(T, 2);
//    }
//    p = 997.18 + 3.1439e-3*T - 3.7574e-3*pow(T, 2);
//
//    T = T + 273.15; /* Convert back to Kelvin to find viscosity */
//    /* Calculate the Reynolds and Prandlt numbers */
//    Re = p*v*L/mu(T);
//    Pr = Cp_wat*mu(T)/k;
//    
//    /* Source: An introduction to Heat and Mass Transfer (Middleman) */
//    return (0.35 + 0.56*pow(Re,0.52))*pow(Pr,0.3)*k/L;
//}

