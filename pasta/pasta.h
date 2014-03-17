#ifndef PASTA_H
#define PASTA_H

#include "constants.h"

/* Composition ****************************************************************/
double conc_gas(double, double, double, double, double);
double conc_air(double, double, double, double, double);
double conc_vap(double, double, double, double, double);

double mdb_wat(double, double, double);
double mdb_wat_sat(double, double);
double sat_wat(double, double, double);
double sat_gas(double, double, double);

double molefrac_vap(double);

/* Fluid **********************************************************************/
double visc_wat(double);
double perm_wat(double, double, double);

/* Gas ************************************************************************/
double rho_gas(double, double, double);
double pvap_wat(double);

/* Phase Change ***************************************************************/
double evap(double, double, double, double, double);

/* Thermal ********************************************************************/
double rho_eff(double, double, double, double, double);
double Cp_eff(double, double, double, double, double);
double k_eff(double, double, double, double, double);
double fluidconv(double, double, double, double, double);

#endif

