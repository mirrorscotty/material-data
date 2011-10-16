#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include "datafile.h"
#include "can.h"

/* Data file parsing functions. Most of these are specific to just the can
 * model because of annoying global variables.
 * TODO: Get rid of them!
 */

extern double Mpro, Mfat, Mcar, Mfib, Mash, Mwat, Mice;
extern double AA, EaA, AB, EaB, R;
extern double A, B, C;
extern double Pressure, molar_mass;
extern double Sutherland, Tref, muref;
extern double To, Text_hot, Text_cold;
extern double v, L, t_heat;

/**
 * Initialization function for the library.
 */
int get_vars(char *str)
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

    delete_buffer(buffer);
    //free(vars);

	return 1;
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
	FIND("AA", data, line)
	FIND("EaA", data, line)
	FIND("AB", data, line)
	FIND("EaB", data, line)
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
	STO(AA, data)
	STO(EaA, data)
	STO(AB, data)
	STO(EaB, data)
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
