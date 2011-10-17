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
double MW_wat, MW_pro, MW_fat, MW_car, MW_fib, MW_ash;
double Mwat, Mpro, Mfat, Mcar, Mfib, Mash, Mice;
double Ea, A, To, Tf, L, R;

int main(int argc, char *argv[])
{
    init("freezing_data.dat");
	output_data();
	return 0;
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
	fp = fopen(file, "w");

	min = 100;
	max = 280;
	points = 100;

	double *x, *y;
	x = (double*) malloc(sizeof(double)*points);
	y = (double*) malloc(sizeof(double)*points);

    //MALLOC_CHECK(x || y)

	for (i=1; i < points; i++) {
		x[i] = min+(max-min)*i/points;
		y[i] = Cp(x[i]);
		fprintf(fp, "%f,%f\n", x[i], y[i]);
	}

	fclose(fp);
//	free(x);
//	free(y);

	return 0;
}

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

	FIND("MW_wat", data, line)
	FIND("MW_pro", data, line)
	FIND("MW_fat", data, line)
	FIND("MW_car", data, line)
	FIND("MW_fib", data, line)
	FIND("MW_ash", data, line)
	FIND("Mwat", data, line)
	FIND("Mpro", data, line)
	FIND("Mfat", data, line)
	FIND("Mcar", data, line)
	FIND("Mfib", data, line)
	FIND("Mash", data, line)
	FIND("Mice", data, line)
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
        STO(data, MW_wat)
        STO(data, MW_pro)
        STO(data, MW_fat)
        STO(data, MW_car)
        STO(data, MW_fib)
        STO(data, MW_ash)
        STO(data, Mwat)
        STO(data, Mpro)
        STO(data, Mfat)
        STO(data, Mcar)
        STO(data, Mfib)
        STO(data, Mash)
        STO(data, Mice)
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
    double MW_s = MW_solids();
	return ( x*MW_wat/(x*MW_wat+(1-x-y)*MW_wat+y*MW_s) );
}

/* Calculate the mole fraction of "a" given it's mass fraction and molecular
 * weight. The total number of mole is calculated from the composition data
 * provided in the data file. */
double MoleFrac(double Ma, double MWa)
{
    double total_moles;
    total_moles = (Mwat+Mice)/MW_wat + Mpro/MW_pro + Mfat/MW_fat + Mcar/MW_car + Mfib/MW_fib + Mash/MW_ash;
    return (Ma/MWa)/total_moles;
}

/* Determine the average molecular weight of the solids. */
double MW_solids()
{
    return (MW_pro + MW_fat + MW_car + MW_fib + MW_ash)/5.0;
}

/* Mole fraction of solids */
double X_solids()
{
    return 1-MoleFrac((Mwat+Mice), MW_wat);
}

/**
 * Calculate the mole fraction of ice in food given the temperature and mole
 * fraction of solids
 */
double X_ice(double T)
{
	/* Calculate the initial freezing temperature */
	double Ti, x_w1, Xs;

    Xs = X_solids();
	Ti = pow( (1/Tf - R/L*log(1-Xs)), -1);
	
	if(T>Ti) {
		return 0; /* No ice formed above the freezing point */
	} else {
		x_w1 = exp((1/Tf - 1/T) * L/R);
		return ( 1-x_w1-X_solids() );
	}
}

/**
 * Calculate the heat capacity of the partially frozen food products.
 */
double Cp(double T)
{
	double Mi, Mw, Ms, dMi, dMw, Ti;
    double Xw, Xs;
	double dT = 0.0001; /* Used for calculating derivatives */

    Xw = MoleFrac(Mwat, MW_wat);
    Xs = X_solids();

	/* Calculate the initial freezing temperature */
	Ti = pow( (1/Tf - R/L*log(1-Xs)), -1);
	
	Mi = M_ice(X_ice(T), Xs); /* Mass fraction of ice */
	Mw = M_ice(Xw-X_ice(T), Xs); /* Mass fraction of water */
	Ms = 1-Mi-Mw; /* Mass fraction of solids */

	dMi = (M_ice(X_ice(T+dT), Xs)-Mi)/dT;
	dMw = (M_ice(Xw-X_ice(T+dT), Xs)-Mw)/dT;

	return ( Mw*Cp_water(T) + Ms*Cp_solids(T) + Mi*Cp_ice(T) - L*dMi -
	         (dMw*Cp_water(T) + dMi*Cp_ice(T))*(Ti-T) );
}

double Cp_water(double T)
{
    T = T-273.15;
    if(T >= 0) {
        return 4.1289 + 9.0864e-5*T - 5.4731e-6*pow(T, 2);
    } else {
        return 4.1289 + 5.3062e-3*T - 9.9516e-4*pow(T, 2);
    }
}

double Cp_solids(double T)
{
    T = T-273.15;
    double Cp_pro, Cp_fat, Cp_car, Cp_ash, Cp_fib;
    Cp_pro = 2.0082 + 1.2089e-3*T - 1.3129e-6*pow(T, 2);
    Cp_fat = 1.9842 + 1.4733e-4*T - 4.8008e-6*pow(T, 2);
    Cp_car = 1.5488 + 1.9625e-3*T - 5.9399e-6*pow(T, 2);
    Cp_ash = 1.0926 + 1.8896e-3*T - 3.6817e-6*pow(T, 2);
    Cp_fib = 1.8459 + 1.8306e-3*T - 4.6509e-6*pow(T, 2);

    return Mpro*Cp_pro + Mfat*Cp_fat + Mcar*Cp_car + Mfib*Cp_fib + Mash*Cp_ash;
}

double Cp_ice(double T)
{
    T = T-273.15;
    return 2.0623 + 6.0769e-3*T;
}

/* Calculate the thermal conductivity using the Choi-Okos equations. */
double k(double T)
{
    T = T-273.15;
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

/* Calculate density using the Choi-Okos equations. */
double rho(double T)
{
    T = T-273.15;
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

double p_solids(double T)
{
    T = T-273.15;
    double p_pro, p_fat, p_car, p_fib, p_ash, p_wat, p_ice;

    /* Calculate the densities */
    p_pro = 1.3299e3 - 5.1840e-1*T;
    p_fat = 9.2559e2 - 4.1757e-1*T;
    p_car = 1.5991e3 - 3.1046e-1*T;
    p_fib = 1.3115e3 - 3.6589e-1*T;
    p_ash = 2.4238e3 - 2.8063e-1*T;

    return 1/(Mpro/p_pro + Mfat/p_fat + Mcar/p_car + Mfib/p_fib + Mash/p_ash);
}

double p_water(double T)
{
    T = T-273.15;
    return 997.18 + 3.1439e-3*T - 3.7574e-3*pow(T, 2);
}

double p_ice(double T)
{
    T = T-273.15;
    return 916.89 - 1.3071e-1*T;
}

/**
 * Determine the volume fraction of water
 */
double Xv_water(double T)
{
	double Mi, Mw, Ms;
    
	Mi = M_ice(X_ice(T), X_solids());
	Mw = M_ice(MoleFrac(Mwat, MW_wat)-X_ice(T), X_solids());
	Ms = 1-Mi-Mw;

    /* FixMe! */
	return (Mw/p_water(T)) / (Mi/p_ice(T) + Mw/p_water(T) + Ms/p_solids(T));
}

/**
 * Calculate the rate of reaction factoring in increase in concentration as a
 * result of the ice crystals forming.
 */
double reaction_rate(double T, double c)
{
	return -A*exp(-Ea/(R*T))*Xv_water(To)/Xv_water(T)*c;
}


