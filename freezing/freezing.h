#ifndef FREEZING_H
#define FREEZING_H

#ifdef FZ_DLLEXPORT
#define EXTFZ_API __declspec(dllexport)
#else
#define EXTFZ_API
#endif

/* Define a macro to shorten code later on. This handles making sure that all
 * of the property functions don't have to worry about messing with arrays
 * and stuff
 */
#define LOOP(x) for(i=0;i<blockSize;i++){outReal[i]=x(inReal[0][i]);}
/* Do the same thing, only for functions that take two arguments */
#define LOOP2(x) for(i=0;i<blockSize;i++){outReal[i]=x(inReal[0][i],inReal[1][i]);}

/* Function prototypes */
double alpha(double);
double Cp(double); /* Verified */
double k(double); /* Verified */
double reaction_rate(double, double); /* Verified */
double rho(double); /* Verified */
double M_ice(double, double); /* Verified implicitely */
double X_ice(double); /* Verified */
double Xv_water(double); /* Verified */

/* New Functions! */
double MoleFrac(double, double);
double MW_solids();
double X_solids();
double Cp_solids(double);
double Cp_water(double);
double Cp_ice(double);
double p_solids(double);
double p_water(double);
double p_ice(double);

EXTFZ_API int init(const char*);
EXTFZ_API const char * getLastError();
EXTFZ_API int eval(const char*, int, const double**, const double**, int, double*, double*);

#endif
