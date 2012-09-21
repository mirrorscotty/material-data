#ifndef CAN_H
#define CAN_H

#ifdef __cplusplus
extern "C" {
#endif

#define MALLOC_CHECK(x) if(x == NULL){report_error("Failed to allocate memory for #x"); exit(1);}

#define alpha(T) (k((T))/(rho((T))*Cp((T))))

/* Function prototypes */
double Cp(double);
double k(double);
double reaction_rate1(double, double);
double reaction_rate2(double, double);
double T_init(double);
double T_ext(double);
//double h(double);
double rho(double); 
double mu(double);

void initialize_variables();
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
double reaction_rate(double, double);
double CpFz(double);
double Cp_water(double);
double Cp_solids(double);
double Cp_ice(double);


#ifdef __cplusplus
}
#endif

#endif
