#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "freezing.h"

#define MALLOC_CHECK(x) if(!(x)){fprintf(stderr, "Failed to allocate memory"); exit(0);}

/* TODO:
 * Add functions for convective heat transfer coefficient
 * Change composition to mass fraction instead of mole fraction
 * Plug memory leaks
 * Double check functions to make sure they're accurate
 */

typedef struct {
    char *name;
    double value;
} variable;

/* Define a variable to hold any error messages returned by the function */
static const char *error = NULL;

/* Annoying global variable definitions. For now, all the model parameters are
 * hard coded into the dll using these variables. This should be changed to
 * allow them to be read in from a data file.
 */
double MW_w, MW_s, Xw, Xs, Cp_solids, Cp_water, Cp_ice, p_solids, p_water;
double p_ice, k_solids, k_water, k_ice, Ea, A, To, Tf, L, R;

/*
int main(int argc, char *argv[])
{
    init("freezing_data.dat");
	output_data();
	return 0;
}
*/

/* Test function to spit out a table of data with the x values in one column
 * and the results of a function the other.
 */
/*
int output_data()
{
	double min, max;
	int points, i;
	const char *file = "out.csv";

	FILE *fp;
	fp = fopen(file, "w");

	min = 100;
	max = 280;
	points = 100;

	double *x, *y;
	x = (double*) malloc(sizeof(double)*points);
	y = (double*) malloc(sizeof(double)*points);

    //MALLOC_CHECK(x || y)

	for (i=1; i <= points; i++) {
		x[i] = min+(max-min)*i/points;
		y[i] = Cp(x[i]);
		fprintf(fp, "%f,%f\n", x[i], y[i]);
	}

	fclose(fp);
//	free(x);
//	free(y);

	return 0;
}
*/

char** read_datafile(char *filename)
{
    FILE *fp;
    int i;
    int j;
    char **buffer;

    fp = fopen(filename, "r");

    /* Only read the first 200 lines */
    buffer = (char **) malloc(sizeof(char*)*200);
    MALLOC_CHECK(buffer)
    for(i=0; i < 200; i++) {
        /* Only read the first 80 characters of each line */
        buffer[i] = (char*) malloc(sizeof(char)*80);
        MALLOC_CHECK(buffer[i])
    }

    /* Read in the file */
    for(i=0; (i<200) || (!feof(fp)); i++) {
        for(j = 0; (j<80) || (!feof(fp)); j++) {
            fscanf(fp, "%c", &buffer[i][j]);
            if(buffer[i][j] == '\0') {
                break;
            }
	    if(buffer[i][j] == '\n') {
	    	break;
	    }
        }
    }
    
    return buffer;
}

/* Deallocate the memory for the buffer created with "read_datafile" */
int delete_buffer(char** buffer)
{
    int i, j;
    for(i=0; i < 200; i++) {
        free(buffer[i]);
    }
    free(buffer);

    return 0;
}

/* Parse the data file contents and return an array of variable names and
 * values. 
 */
#define FIND(X,Y,STR) if(strstr(STR, X)) {strcpy(Y.name, X);}
variable read_line(char* line)
{
	char *value;

	variable data;
	data.name = (char*) malloc(sizeof(char)*20);

	FIND("MW_w", data, line)
	FIND("MW_s", data, line)
	FIND("Xw", data, line)
	FIND("Xs", data, line)
	FIND("Cp_solids", data, line)
	FIND("Cp_water", data, line)
	FIND("Cp_ice", data, line)
	FIND("p_water", data, line)
	FIND("p_solids", data, line)
	FIND("p_ice", data, line)
	FIND("k_solids", data, line)
	FIND("k_water", data, line)
	FIND("k_ice", data, line)
	FIND("Ea", data, line)
	FIND("A", data, line)
	FIND("To", data, line)
	FIND("Tf", data, line)
	FIND("L", data, line)
	FIND("R", data, line)

	if(strcmp(data.name, "")) {
		if(value = strpbrk(line, "0123456789.-")) {
			data.value = atof(value);
		}
	} else {
		strcpy(data.name,"NULL");
		data.value = 3.141592654;
	}
	if(strcmp(data.name, "NULL") != 0) {
		printf("%s --> %f\n", data.name, data.value);
	}
	
	return data;
}

/* Rather than delete the last part of the string, just replace the first '#'
 * with a null character.
 */
char* remove_comments(char* line)
{
	char* comment;
	comment = strchr(line, '#');
	
	if(comment) {
	//	printf("%s", comment);
		*comment = '\0';
	}

	return line;
}

/* Lets play abuse the preprocessor! */
#define STO(VAR, NAME) if(strcmp(VAR.name, #NAME) == 0) { NAME = VAR.value; }

/* Store the data that has been parsed into the annoyingly ugly global varibles */
/* Also, this function is a terrible hack. */
int store_data(variable data)
{
        STO(data, MW_w)
        STO(data, MW_s)
        STO(data, Xw)
        STO(data, Xs)
        STO(data, Cp_solids)
        STO(data, Cp_water)
        STO(data, Cp_ice)
        STO(data, p_solids)
        STO(data, p_water)
        STO(data, p_ice)
        STO(data, k_solids)
        STO(data, k_water)
        STO(data, Ea)
        STO(data, A)
        STO(data, To)
        STO(data, Tf)
        STO(data, L)
        STO(data, R)

    return 0;
}

#define PRNT(X) printf("Value of %s: %f\n", #X, X);
int print_global_vars()
{
	PRNT(MW_w)
        PRNT(MW_s)
        PRNT(Xw)
        PRNT(Xs)
        PRNT(Cp_solids)
        PRNT(Cp_water)
        PRNT(Cp_ice)
        PRNT(p_solids)
        PRNT(p_water)
        PRNT(p_ice)
        PRNT(k_solids)
        PRNT(k_water)
        PRNT(Ea)
        PRNT(A)
        PRNT(To)
        PRNT(Tf)
        PRNT(L)
        PRNT(R)

	return 0;
}

/* The following three functions are required in order to interface correctly
 * with Comsol
 */
/**
 * Initialization function for the library.
 */
EXTFZ_API int init(const char *str)
{
    //Commented because it's not ready!
    char** buffer;
    //variable vars[20], tmp;
    int i, j;
    j = 0;

    buffer = read_datafile(str);

    for(i=0; i < 200; i++) {
        buffer[i] = remove_comments(buffer[i]);
        store_data(read_line(buffer[i]));
        //if(!strcmp("NULL", tmp.name)) {
        //    vars[j] = tmp;
        //    j++;
        //}
    }


    //vars = parse_datafile(buffer);

	print_global_vars();

    //delete_buffer(buffer);
    //free(vars);


	return 1;
}

/**
 * Error reporting function. Returns the last error message set.
 */
EXTFZ_API const char * getLastError()
{
	return error;
}

/**
 * Function to that actually does the work for the library.
 */
EXTFZ_API int eval(const char *func,
		   int nArgs,
		   const double **inReal,
		   const double **inImag,
		   int blockSize,
		   double *outReal,
		   double *outImag)
{
	int i, j;
	if(strcmp(func, "Cp") == 0) {
		LOOP(Cp)
	} else if(strcmp(func, "rho") == 0) {
		LOOP(rho)
	} else if(strcmp(func, "reaction_rate") == 0) {
		LOOP2(reaction_rate)
	} else if(strcmp(func, "k") == 0) {
		LOOP(k)
	} else {
		error = "Cannot find function";
		return 0;
	}
    return 1;
}

/* Functions to actually calculate stuff. */

/**
 * Calculate the mass fraction of ice given the mole fraction of ice and solids
 * x -> ice/water
 * y -> solids
 */
double M_ice(double x, double y)
{
	return ( x*MW_w/(x*MW_w+(1-x-y)*MW_w+y*MW_s) );
}

/**
 * Calculate the mole fraction of ice in food given the temperature and mole
 * fraction of solids
 */
double X_ice(double T)
{
	/* Calculate the initial freezing temperature */
	double Ti, x_w1;
	Ti = pow( (1/Tf - R/L*log(1-Xs)), -1);
	
	if(T>Ti) {
		return 0; /* No ice formed above the freezing point */
	} else {
		x_w1 = exp((1/Tf - 1/T) * L/R);
		return ( 1-x_w1-Xs );
	}
}

/**
 * Calculate the heat capacity of the partially frozen food products.
 */
double Cp(double T)
{
	double Mi, Mw, Ms, dMi, dMw, Ti;
	double dT = 0.0001; /* Used for calculating derivatives */

	/* Calculate the initial freezing temperature */
	Ti = pow( (1/Tf - R/L*log(1-Xs)), -1);
	
	Mi = M_ice(X_ice(T), Xs); /* Mass fraction of ice */
	Mw = M_ice(Xw-X_ice(T), Xs); /* Mass fraction of water */
	Ms = 1-Mi-Mw; /* Mass fraction of solids */

	dMi = (M_ice(X_ice(T+dT), Xs)-Mi)/dT;
	dMw = (M_ice(Xw-X_ice(T+dT), Xs)-Mw)/dT;

	return ( Mw*Cp_water + Ms*Cp_solids + Mi*Cp_ice - L*dMi -
	         (dMw*Cp_water + dMi*Cp_ice)*(Ti-T) );
}

/**
 * Calculate thermal conductivity as a function of temperature
 */
double k(double T)
{
	double Mi, Mw, Ms, Xvi, Xvs, Xvw;
	Mi = M_ice(X_ice(T), Xs);
	Mw = M_ice(Xw-X_ice(T), Xs);
	Ms = 1-Mi-Mw;

	Xvi = (Mi/p_ice) / (Mi/p_ice + Mw/p_water + Ms/p_solids);
	Xvw = (Mw/p_water) / (Mi/p_ice + Mw/p_water + Ms/p_solids);
	Xvs = (Ms/p_solids) / (Mi/p_ice + Mw/p_water + Ms/p_solids);

	return k_ice*Xvi + k_water*Xvw + k_solids*Xvs;
}

/**
 * Calculate the density of the slab as it freezes as a function of temperature
 * using the Choi-Okos equations.
 */
double rho(double T)
{
	double Mi, Mw, Ms;
	Mi = M_ice(X_ice(T), Xs);
	Mw = M_ice(Xw-X_ice(T), Xs);
	Ms = 1-Mi-Mw;
	
	return 1/(Mi/p_ice + Mw/p_water + Ms/p_solids);
}

/**
 * Determine the volume fraction of water
 */
double Xv_water(double T)
{
	double Mi, Mw, Ms;
	Mi = M_ice(X_ice(T), Xs);
	Mw = M_ice(Xw-X_ice(T), Xs);
	Ms = 1-Mi-Mw;
	
	return (Mw/p_water) / (Mi/p_ice + Mw/p_water + Ms/p_solids);
}

/**
 * Calculate the rate of reaction factoring in increase in concentration as a
 * result of the ice crystals forming.
 */
double reaction_rate(double T, double c)
{
	return -A*exp(-Ea/(R*T))*Xv_water(To)/Xv_water(T)*c;
}


