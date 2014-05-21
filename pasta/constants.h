/**
 * @file constants.h
 * Set of constants for calculating pasta properties
 */

/**
 * Ideal gas constant [J/(K mol)]
 */
#define GASCONST 8.3144621
/**
 * Molar mass of water [g/mol]
 */
#define MWWAT 18.01528
/**
 * Molar mass of air [g/mol]
 */
#define MWAIR 28.97

/**
 * Composition for pure water, for use in CreateChoiOkos. Components are Mpro,
 * Mfat, Mcar, Mfib, Mash, Mwat, Mice (all mass fractions)
 *
 * @see CreateChoiOkos
 */
#define WATERCOMP 0, 0, 0, 0, 0, 1, 0 /* Composition of water */
/**
 * Composition of pasta From Cummings et al. 1993, no water. For use in
 * CreateChoiOkos Components are Mpro, Mfat, Mcar, Mfib, Mash, Mwat, Mice (all
 * mass fractions)
 *
 * @see CreateChoiOkos
 */
#define PASTACOMP .0166, .0239, .784, .0179, .00824, 0, 0

/* Heat capacities of random stuff */
/**
 * Vapor heat capacity at 100C [J/(kg K)]
 */ 
#define CPVAPOR 2080
/**
 * Heat capacity of air at 23C, 760 mmHg [J/(kg K)]
 */
#define CPAIR 1012

/* Thermal conductivities */
/**
 * Thermal conductivity for air at 300K, 1 bar [W/(m K)]
 */
#define KAIR 0.0262
/**
 * Vapor thermal conductivity at 293 K, 101.3 kPa [W/(m K)]
 */
#define KVAPOR 0.016

/**
 * Intrinsic permeability of water [m^2]  Source: Zhu et al. 2011
 */
#define PERMWAT 1e-18

/**
 * Intrinsic permeability gas [m^2]
 */
#define PERMGAS 0

/* Random extra values */
/**
 * Porosity [-] (From Xiong et al. 1991)
 */
#define POROSITY 0.0612
/**
 * Intrinsic Permeability From Zhu et al. 2011
 */
#define KWINTR 1e-18
/**
 * Irreducible water saturation Note: This value is made up
 */
#define SR .02
