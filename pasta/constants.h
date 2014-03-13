#define GASCONST 8.3144621 /* Ideal gas constant [J/(K mol)] */
#define MWWAT 18.01528 /* Molar mass of water [g/mol] */
#define MWAIR 28.97 /* Molar mass of air [g/mol] */

/* Mpro, Mfat, Mcar, Mfib, Mash, Mwat, Mice */
#define WATERCOMP 0, 0, 0, 0, 0, 1, 0 /* Composition of water */
#define PASTACOMP 0, 0, 0, 1, 0, 0, 0 /* Composition of pasta */

/* Heat capacities of random stuff */
#define CPVAPOR 2080 /* Vapor heat capacity at 100C [J/(kg K)] */
#define CPAIR 1012 /* Heat capacity of air at 23C, 760 mmHg [J/(kg K)] */

/* Thermal conductivities */
#define KAIR 0.0262 /* Thermal conductivity for air at 300K, 1 bar [W/(m K)] */
#define KVAPOR 0.016 /* 293 K, 101.3 kPa [W/(m K)] */

/* Intrinsic permeability */
#define PERMWAT 1e-18 /* From Zhu et al. 2011 */

/* Porosity */
#define POROSITY 0.2 /* Selected arbitrarily */

