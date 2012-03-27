#ifndef CAN_H
#define CAN_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef CAN_DLLEXPORT
#define EXTCAN_API __declspec(dllexport)
#else
#define EXTCAN_API
#endif

/* Define a macro to shorten code later on. This handles making sure that all
 * of the property functions don't have to worry about messing with arrays
 * and stuff
 */
#define LOOP(x) for(i=0;i<blockSize;i++){outReal[i]=x(inReal[0][i]);}
/* Do the same thing, only for functions that take two arguments */
#define LOOP2(x) for(i=0;i<blockSize;i++){outReal[i]=x(inReal[0][i],inReal[1][i]);}
#define MALLOC_CHECK(x) if(x == NULL){report_error("Failed to allocate memory for #x"); exit(1);}


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

EXTCAN_API int init(const char*);
EXTCAN_API const char * getLastError();
EXTCAN_API int eval(const char*, int, const double**, const double**, int, double*, double*);

#ifdef __cplusplus
}
#endif

#endif
