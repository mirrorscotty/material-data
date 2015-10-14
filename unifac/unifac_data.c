#include "unifac.h"
#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix/2dmatrix/xstrtok.h"

#define UNIFAC_MAXROWS 100
#define UNIFAC_LINELENGTH 100

unifac_data* UnifacCreateData()
{
    unifac_data *d;
    d = (unifac_data*) calloc(sizeof(unifac_data), 1);
    return d;
}

void UnifacPrintRow(unifac_row d)
{
    printf("%d\t%d\t%s\t%s\t\t%.4f\t%.4f\n",
            d.id, d.groupid, d.group, d.subgroup, d.Rk, d.Qk);
    return;
}

void UnifacPrintTable(unifac_data *d)
{
    int i;
    printf("ID\tGID\tGroup\tSubgroup\tRk\tQk\n");
    printf("------------------------------------------------------\n");
    for(i=0; i<d->ngroups; i++)
        UnifacPrintRow(d->rows[i]);

    puts("");
    mtxprnt(d->interactions);
    return;
}

void UnifacLoadData(unifac_data *d, char* filename)
{
    int maxlines = UNIFAC_MAXROWS,
        maxchars = UNIFAC_LINELENGTH,
        i, j,
        ncols = 5, /* Assume there's at least one column */
        nrows = 0,
        row0 = 1; /* Include a single line header */
    const char *delim = ",";
    char **buffer;
    char *number;
    FILE *fp; /* TODO: Actually open/read the file */

    /* Make a buffer to store all of the characters read from the file */
    buffer = (char**) calloc(sizeof(char*), maxlines);
    for(i=0; i<maxlines; i++)
        buffer[i] = (char*) calloc(sizeof(char), maxchars);

    /* Load in all the lines (up to the specified limits) */
    fp = fopen(filename, "r");
    for(i=0; i<maxlines; i++) {
        if(fgets(buffer[i], maxchars, fp) == NULL) {
            nrows = i;
            break;
        }
    }
    if(i==0)
        fprintf(stderr, "Maximum number of lines exceeded. Current limit is %d lines.\n", maxlines);
    fclose(fp);

    /* Check the first line to see how many values there are. */
    i=0;
    while(buffer[row0][i]) {
        if(buffer[row0][i] == delim[0])
            ncols++;
        i++;
    }

    nrows = nrows - row0;

    /* Make a matrix that's hopefully the right size */
    d->rows = (unifac_row*) calloc(sizeof(unifac_row), nrows);
    d->ngroups = nrows;

    /* Start putting values into it */
    for(i=0; i<nrows; i++) {
        /* Set the subgroup ID */
        d->rows[i].id = i;

        /* Get the group ID and store it. */
        number = xstrtok(buffer[i+row0], delim);
        if(number[0] == '\0')
            d->rows[i].groupid = -1;
        else
            d->rows[i].groupid = atof(number);

        /* Group name */
        number = xstrtok(NULL, delim);
        strncpy(d->rows[i].group, number, 20);
        number = xstrtok(NULL, delim);
        strncpy(d->rows[i].subgroup, number, 20);

        number = xstrtok(NULL, delim);
        if(number[0] == '\0')
            d->rows[i].Rk = -1;
        else
            d->rows[i].Rk = atof(number);

        number = xstrtok(NULL, delim);
        if(number[0] == '\0')
            d->rows[i].Qk = -1;
        else
            d->rows[i].Qk = atof(number);
    }

    for(i=0; i<maxlines; i++)
        free(buffer[i]);
    free(buffer);
    //free(number);

    return;
}

void UnifacLoadInteractions(unifac_data *d, char *filename)
{
    d->interactions = mtxloadcsv(filename, 1);
    return;
}

