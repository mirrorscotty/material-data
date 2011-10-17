#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "freezing.h"

#define MALLOC_CHECK(x) if(!(x)){fprintf(stderr, "Failed to allocate memory"); exit(0);}

typedef struct {
    char name[20];
    double value;
} variable;

/* Define a variable to hold any error messages returned by the function */
static const char *error = NULL;

/* Annoying global variable definitions. For now, all the model parameters are
 * hard coded into the dll using these variables. This should be changed to
 * allow them to be read in from a data file.
 */

double MW_w = 18; // Molecular weight of water
double MW_s = 2.6148e5; // Average molecular weight of solids

double Xw = .9966; // Initial mole fraction of water
double Xs = .0034; // Mole fraction of solids

// Heat capacities in J/(mol K)
double Cp_solids = 50;
double Cp_water = 75.3;
double Cp_ice = 60;

// Densities in kg/m^3
double p_solids = 1200;
double p_water = 1000;
double p_ice = 916.9;

// Heat capacities in W/(m K)
double k_solids = 3;
double k_water = 0.596;
double k_ice = 2.1;

// Rate of reaction parameters
double Ea = 10000;
double A = 1;

double T0 = 280; // Initial Temperature
double Tf = 273.15; // Freezing point of pure water

double L = 6010; // J/mol Latent heat of fusion of water
double R = 8.3145; // J/(mol K)




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

	for (i=1; i <= points; i++) {
		x[i] = min+(max-min)*i/points;
		y[i] = Cp(x[i]);
		fprintf(fp, "%f,%f\n", x[i], y[i]);
	}

	fclose(fp);
	free(x);
	free(y);

	return 0;
}


char** read_datafile(char *filename)
{
    FILE *fp;
    fp = fopen(filename, "r");
    int i, j;

    char **buffer;
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
}

/* Parse the data file contents and return an array of variable names and
 * values. 
 */
variable* parse_datafile(char **buffer)
{
    /* Allocate space for 19 variables total */
    variable* lines;
    lines = (variable*) malloc(sizeof(variable)*19);
    MALLOC_CHECK(lines)
    
    int i, j;
    int nameread = 0;
    int psn = 0;
    int var = 0; /* Keep track of the current variable */
    char* tmp;
    tmp = (char*) malloc(sizeof(char)*20);
    MALLOC_CHECK(tmp)

    for(i=0; i<200; i++) {
        psn = 0;
        for(j=0; j<80; j++) {
            /* Ignore comments and the end of the line */
            if((buffer[i][j] == '#') ||
               (buffer[i][j] == '\0') ||
               (buffer[i][j] == '\n')) {
               tmp[psn] = '\0'; /* Null terminate the tmp buffer so nothing */
                                /* crazy happens. */
               break;
            }
            /* If we haven't read in the variable's name yet, do that. */
            if(!nameread) {
                /* If the character is an equals sign, then we've got the name */
                if (buffer[i][j] != '=') {
                    /* Don't read spaces */
                    if (buffer[i][j] != ' ') {
                        tmp[psn] = buffer[i][j];
                        psn++;
                    }
                } else {
                    tmp[psn+1] = '\0'; /* Null terminate the string */
                    nameread = 1; /* Name: read */
                    psn = 0; /* Reset the position in the tmp variable to 0 */
                    strcpy(lines[var].name, tmp);
                }
            } else { /* Read in the number now */
                /* Still not reading spaces */
                if (buffer[i][j] != ' ') {
                    tmp[psn] = buffer[i][j];
                    psn++;
                }
            }
        }
        /* If a variable name was read, then assume that the value was also
         * read. Store the value and go to the next variable.
         */
        if (nameread) {
            /* Convert the value to a double and store it */
            lines[var].value = atof(tmp); 
            var++;
        }
    }

    //free(tmp);

    return lines;
}

/* Lets play abuse the preprocessor! */
#define STO(VAR, NAME) if(strcmp(VAR.name, "NAME")) { NAME = VAR.value; }

/* Store the data that has been parsed into the annoyingly ugly global varibles */
/* Also, this function is a terrible hack. */
int store_data(variable* data)
{
    int i;
    for(i=0; i<19; i++) {
        STO(data[i], MW_w)
        STO(data[i], MW_s)
        STO(data[i], Xw)
        STO(data[i], Xs)
        STO(data[i], Cp_solids)
        STO(data[i], Cp_water)
        STO(data[i], Cp_ice)
        STO(data[i], p_solids)
        STO(data[i], p_water)
        STO(data[i], p_ice)
        STO(data[i], k_solids)
        STO(data[i], k_water)
        STO(data[i], Ea)
        STO(data[i], A)
        STO(data[i], T0)
        STO(data[i], Tf)
        STO(data[i], L)
        STO(data[i], R)
    }
}

/**
 * Initialization function for the library.
 */
EXTFZ_API int init(const char *str)
{
    //Commented because it's not ready!
    char** buffer;
    variable* vars;

    buffer = read_datafile(str);
    vars = parse_datafile(buffer);
    store_data(vars);

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
	return -A*exp(-Ea/(R*T))*Xv_water(T0)/Xv_water(T)*c;
}


