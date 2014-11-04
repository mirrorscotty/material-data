#ifndef DRYING_H
#define DRYING_H

#ifdef CAN_DLLEXPORT
#define EXTCAN_API __declspec(dllexport)
#else
#define EXTCAN_API
#endif

#include "regress.h"

/* Define a macro to shorten code later on. This handles making sure that all
 * of the property functions don't have to worry about messing with arrays
 * and stuff
 */
#define LOOP(x) for(i=0;i<blockSize;i++){outReal[i]=x(inReal[0][i]);}
/* Do the same thing, only for functions that take two arguments */
#define LOOP2(x) for(i=0;i<blockSize;i++){outReal[i]=x(inReal[0][i],inReal[1][i]);}

/* Function prototypes */
void initialize_variables();
int output_data();
matrix* create_D_data();

EXTCAN_API int init(const char*);
EXTCAN_API const char * getLastError();
EXTCAN_API int eval(const char*, int, const double**, const double**, int, double*, double*);

#endif
