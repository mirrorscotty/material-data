#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "drying.h"
#include "drying_D.h"
#include "regress.h"

#ifdef __linux__ /* Not using windows */
#include<errno.h> /* Include the standard error reporting library */
#endif

#define MALLOC_CHECK_ 1

#define MALLOC_CHECK(x) if(x == NULL){report_error("Failed to allocate memory."); exit(1);}

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

double Mpro, Mfat, Mcar, Mfib, Mash, Mwat, Mice;
double A1, Ea1, A2, Ea2, R;
double A, B, C;
double Pressure, molar_mass;
double Sutherland, Tref, muref;
double To, Text_hot, Text_cold;
double v, L, t_heat;

/* For testing only */
int main(int argc, char *argv[])
{
    int i=0;
    set_default_globals();
    DestroyMatrix(create_D_data());
}

/* Set all global variables to an initial value of zero in case something goes
 * horrible, horribly wrong and something tries to read an uninitialized value.
 */
 /* Needs to be revised/moved to drying_D.c */
void initialize_variables()
{
    Mpro = 0;
    Mfat = 0;
    Mcar = 0;
    Mfib = 0;
    Mash = 0;
    Mwat = 0;
    Mice = 0;
    A1 = 0;
    Ea1 = 0;
    A2 = 0;
    Ea2 = 0;
    R = 0;
    A = 0;
    B = 0;
    C = 0;
    Pressure = 0;
    molar_mass = 0;
    Sutherland = 0;
    Tref = 0;
    muref = 0;
    To = 0;
    Text_hot = 0;
    Text_cold = 0;
    v = 0;
    L = 0;
    t_heat = 0;
}

int report_error(const char *str)
{
    #ifdef __linux__
    fprintf(stderr, "%s: %s\n", str, strerror(errno));
    return errno;
    #else
    fprintf(stderr, "Fatal error: %s\n", str);
    return -1
    #endif
}

/* Test function to spit out a table of data with the x values in one column
 * and the results of a function the other.
 */
/* Needs to be revised */
int output_data()
{
	double min, max;
	int points, i;
	const char *file = "out.csv";

	FILE *fp;
	double *x, *y;

    min = 0;
    max = 0;
    points = 0;
    i = 0;
    fp = NULL;
    x = NULL;
    y = NULL;

	fp = fopen(file, "w");
    if(!fp) {
        report_error("Failed to open file for writing.");
    }

	min = 300;
	max = 480;
	points = 100;

	x = (double*) malloc(sizeof(double)*points);
	y = (double*) malloc(sizeof(double)*points);

    if(x == NULL) {
        report_error("Cannot allocate memory for x");
        exit(1);
    }
    if(y == NULL) {
        report_error("Cannot allocate memory for y");
        exit(1);
    }

    MALLOC_CHECK(x || y)

	for (i=1; i <= points; i++) {
		x[i] = min+(max-min)*i/points;
		//y[i] = Cp(x[i]);
		fprintf(fp, "%f,%f\n", x[i], y[i]);
	}

    if(fp) {
        if(fclose(fp) != 0) {
           report_error("Failed to close file.");
           exit(1);           
        }
    }

	free(x);
	free(y);

	return 0;
}


char** read_datafile(char *filename)
{
    FILE *fp;
    int i, j;
    char **buffer;

    fp = NULL;
    i = 0;
    j = 0;
    buffer = NULL;

    fp = fopen(filename, "r");
    if(!fp) {
       report_error("Cannot open data file for reading.");
       exit(1);
    }
       
    /* Only read the first 200 lines */
    buffer = (char **) malloc(sizeof(char*)*200);
    MALLOC_CHECK(buffer)
    for(i=0; i < 200; i++) {
        /* Only read the first 80 characters of each line */
        buffer[i] = (char*) malloc(sizeof(char)*80);
        memset(buffer[i], '\0', 80);
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
    
    if(fp) {
        if(fclose(fp) != 0) {
           report_error("Failed to close data file.");
           exit(1);           
        }
    }

    return buffer;
}

/* Deallocate the memory for the buffer created with "read_datafile" */
int delete_buffer(char** buffer)
{
    int i = 0;
    for(i=0; i < 200; i++) {
        free(buffer[i]);
    }
    free(buffer);

    return 0;
}

/* Parse the data file contents and return an array of variable names and
 * values. The macro is horrible and uses the worst bit of code imaginable:
 * the goto statment.
 */
#define FIND(X,Y,STR) if(strstr(STR, X)) {strcpy(Y.name, X); goto end;}
variable read_line(char* line)
{
	char *value;
	variable data;

    value = NULL;
    data.name = NULL;
    data.value = 0;
	data.name = (char*) malloc(sizeof(char)*20);
    memset(data.name, '\0', 20);

    if(data.name == NULL) {
        report_error("Null pointer error.");
    }

	FIND("Mpro", data, line)
	FIND("Mfat", data, line)
	FIND("Mcar", data, line)
	FIND("Mfib", data, line)
	FIND("Mash", data, line)
	FIND("Mwat", data, line)
	FIND("Mice", data, line)
	FIND("A1", data, line)
	FIND("Ea1", data, line)
	FIND("A2", data, line)
	FIND("Ea2", data, line)
	FIND("R", data, line)
	FIND("A", data, line)
	FIND("B", data, line)
	FIND("C", data, line)
	FIND("Pressure", data, line)
	FIND("molar_mass", data, line)
	FIND("Sutherland", data, line)
	FIND("Tref", data, line)
	FIND("muref", data, line)
	FIND("To", data, line)
	FIND("Text_hot", data, line)
	FIND("Text_cold", data, line)
	FIND("v", data, line)
	FIND("t_heat", data, line)
	FIND("L", data, line)

    end:

	if(strcmp(data.name, "")) {
		if((value = strpbrk(line, "0123456789.-"))) {
			data.value = atof(value);
		}
	} else {
		strcpy(data.name, "NULL");
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
    comment = NULL;
	comment = strchr(line, '#');
	
	if(comment) {
	//	printf("%s", comment);
		*comment = '\0';
	}

	return line;
}

/* Lets play abuse the preprocessor! */
#define STO(NAME, VAR) if(strcmp(VAR.name, #NAME) == 0) { NAME = VAR.value; }

/* Store the data that has been parsed into the annoyingly ugly global varibles */
/* Also, this function is a terrible hack. */
int store_data(variable data)
{
 	STO(Mpro, data)
	STO(Mfat, data)
	STO(Mcar, data)
	STO(Mfib, data)
	STO(Mash, data)
	STO(Mwat, data)
	STO(Mice, data)
	STO(A1, data)
	STO(Ea1, data)
	STO(A2, data)
	STO(Ea2, data)
	STO(R, data)
	STO(A, data)
	STO(B, data)
	STO(C, data)
	STO(Pressure, data)
	STO(molar_mass, data)
	STO(Sutherland, data)
	STO(Tref, data)
	STO(muref, data)
	STO(To, data)
	STO(Text_hot, data)
	STO(Text_cold, data)
	STO(v, data)
	STO(t_heat, data)
	STO(L, data)

    /* Free the memory allocated to store the variable's name. */
    free(data.name);

    return 0;
}

#define PRNT(X) printf("Value of %s: %f\n", #X, X);
int print_global_vars()
{
    /* TODO: Fix this function */
	PRNT(R)

	return 0;
}

/* The following three functions are required in order to interface correctly
 * with Comsol
 */
/**
 * Initialization function for the library.
 */
EXTCAN_API int init(const char *str)
{
    //Commented because it's not ready!
    char** buffer;
    //variable vars[20], tmp;
    int i, j;
    i = 0;
    j = 0;
    buffer = NULL;

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
EXTCAN_API const char * getLastError()
{
	return error;
}

/**
 * Function to that actually does the work for the library.
 */
EXTCAN_API int eval(const char *func,
		   int nArgs,
		   const double **inReal,
		   const double **inImag,
		   int blockSize,
		   double *outReal,
		   double *outImag)
{
	int i;
	if(strcmp(func, "Cp") == 0) {
        //LOOP(mu)
    } else {
		error = "Cannot find function";
		return 0;
	}
    return 1;
}

/* Functions to actually calculate stuff. */
matrix* create_D_data()
{
    int x_pts, y_pts, i, j;
    int por = .2;
    matrix *X, *Y, *x, *y, *z, *beta;

    x_pts = 10;
    y_pts = 10;
    x = linspace(280, 350, x_pts);
    y = linspace(.05, .9, y_pts);
    z = CreateMatrix(x_pts, y_pts);
    X = CreateMatrix(x_pts*y_pts, 4);
    Y = CreateMatrix(x_pts*y_pts, 1);

    for(i=0; i<x_pts; i++) {
        for(j=0; j<x_pts; j++) {
            setval(z, drying_D(val(x, 0, i), val(y, 0, j), por), i, j);
        }
    }

    for(i=0; i<x_pts; i++) {
        for(j=0; j<y_pts; j++) {
            setval(X, 1, i*y_pts+j, 0);
            setval(X, val(x, 0, i), i*y_pts+j, 1);
            setval(X, val(y, 0, j), i*y_pts+j, 2);
            setval(X, val(x, 0, i)*val(y, 0, j), i*y_pts+j, 3);
            setval(Y, val(z, i, j), i*y_pts+j, 0);
        }
    }

    beta = regress(Y, X);
    //mtxprnt(beta);
    
    DestroyMatrix(x);
    DestroyMatrix(y);
    DestroyMatrix(z);
    DestroyMatrix(X);
    DestroyMatrix(Y);

    return beta;
}

double CalcPoly(matrix* beta, double x, double y)
{
    double a, b, c, d;
    
    a = val(beta, 0, 0);
    b = val(beta, 1, 0);
    c = val(beta, 2, 0);
    d = val(beta, 3, 0);

    return a+x*b+c*y+d*x*y;
}

