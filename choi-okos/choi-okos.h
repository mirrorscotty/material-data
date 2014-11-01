/**
 * @file choi-okos.h
 * Define all the data structures used in the Choi-Okos equations.
 */

#ifndef CAN_H
#define CAN_H

#ifdef __cplusplus
extern "C" {
#endif

#define MALLOC_CHECK(x) if(x == NULL){report_error("Failed to allocate memory for #x"); exit(1);}

#define alpha(comp, T) (k((comp), (T))/(rho((comp), (T))*Cp((comp), (T))))

/**
 * Structure to hold all of the mass fractions and other composition data for
 * food samples. Used for calculation of heat transfer coefficient, density, and
 * heat capacity.
 */
typedef struct {
    double R; /**< Gas constant */

    /* Mass fractions of components */
    double Mpro; /**< Protein */
    double Mfat; /**< Fat */
    double Mcar; /**< Carbohydrates */
    double Mfib; /**< Fiber */
    double Mash; /**< Ash */
    double Mwat; /**< Water */
    double Mice; /**< Ice */

    /* Molar masses of pure components */
    double MW_pro;
    double MW_fat;
    double MW_car;
    double MW_fib;
    double MW_ash;
    double MW_wat;

    /* Freezing parameters */
    double Hfus; /**< Heat of fusion (for water) */
    double Hvap; /**< Heat of vaporization (for water) */
    double Tf; /**< Normal freezing point (for water) */
} choi_okos;

/* Function prototypes */
double alphaFZ(double);
double Cp(choi_okos*, double);
double k(choi_okos*, double);
double rho(choi_okos*, double); 
double mu(choi_okos*, double);

choi_okos *CreateChoiOkos(double, double, double, double, double, double, double);
choi_okos* AddDryBasis(choi_okos*, double);
void DestroyChoiOkos(choi_okos*);
int output_data();
int report_error(const char*);

double M_ice(double, double);
double MoleFrac(double, double);
double MW_solids();
double X_solids();
double X_ice(double);
double p_solids(double);
double p_water(double);
double p_ice(double);
double Xv_water(double);
double Xv_ice(double);
double CpFz(double);
double Cp_water(double);
double Cp_solids(double);
double Cp_ice(double);
double IceMassFrac(double);


#ifdef __cplusplus
}
#endif

#endif

