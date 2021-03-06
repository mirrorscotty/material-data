#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "regress.h"

double square(double x)
{
    return pow(x, 2);
}

matrix* regress(matrix *y, matrix *X)
{
    matrix *Xt, *XtX, *XtXinv, *XtXinvXt, *beta;

    Xt = mtxtrn(X);
    XtX = mtxmul(Xt, X);
    XtXinv = CalcInv(XtX);
    XtXinvXt = mtxmul(XtXinv, Xt);

    beta = mtxmul(XtXinvXt, y);

    DestroyMatrix(Xt);
    DestroyMatrix(XtX);
    DestroyMatrix(XtXinv);
    DestroyMatrix(XtXinvXt);

    return beta;
}

matrix* linspace(double start, double end, int nelem)
{
    matrix* x;
    int i;
    
    x = CreateMatrix(1, nelem);
    for(i=0; i<nelem; i++) {
        setval(x, start + i*(end-start)/(nelem-1), 0, i);
    }

    return x;
}
    
matrix* polyfit(matrix* x, matrix* y, int order)
{
    matrix *Y, *X, *beta;
    int i, j, nelem;

    X = NULL;
    nelem = nCols(x); /* Row matrix */

    X = CreateMatrix(nelem, order+1);

    for(i=0; i<=order; i++) {
        for(j=0; j<nelem; j++) {
            setval(X, pow(val(x, 0, j), i), j, i);
        }
    }
    Y = mtxtrn(y);

    beta = regress(Y, X);
    DestroyMatrix(X);
    DestroyMatrix(Y);

    return beta;
}
 
/* Get the length of a 1D matrix of doubles */
int nCols(matrix *A)
{
    return A->cols;
}

/* Return the number of rows in a two dimensional matrix of doubles */
int nRows(matrix *A)
{
    return A->rows;
}

double val(matrix *A, int row, int col)
{
    if(row >= nRows(A)) {
        fprintf(stderr, "Error: index out of bounds.\n");
        return;
    }
    if(col >= nCols(A)) {
        fprintf(stderr, "Error: Index out of bounds.\n");
        return;
    }
    return A->array[row][col];
}

void setval(matrix *A, double value, int row, int col)
{
    if(row >= nRows(A)) {
        fprintf(stderr, "Error: index out of bounds.\n");
        return;
    }
    if(col >= nCols(A)) {
        fprintf(stderr, "Error: Index out of bounds.\n");
        return;
    }
    A->array[row][col] = value;
}

void mtxprnt(matrix *A)
{
    int i, j;
    
    for(i=0; i<nRows(A); i++) {
        printf("[ ");
        for(j=0; j<nCols(A); j++) {
            printf("%e ", val(A, i, j));
        }
        printf("]\n");
    }
}

void mtxprntfile(matrix *A, char *filename)
{
    int i, j;
    FILE *file;

    file = fopen(filename, "w");
    
    for(i=0; i<nRows(A); i++) {
        for(j=0; j<nCols(A); j++) {
            fprintf(file, "%e,", val(A, i, j));
        }
	fprintf(file, "\n");
    }
    fclose(file);
}

/* Return the transpose of a square matrix */
matrix* mtxtrn(matrix *x)
{
    matrix *xt;
    int rows = nRows(x);
    int cols = nCols(x);
    int i, j;

    xt = NULL;
/*
    if(nRows(x) != nCols(x)) {
        fprintf(stderr, "ERROR!");
        return xt;
    }
*/

    xt = CreateMatrix(cols, rows);

    for(i=0; i<cols; i++) {
        for(j=0; j<rows; j++) {
            setval(xt, val(x, j, i), i, j);
        }
    }

    return xt;
}

/* Multiply matricies using nifty index notation */
matrix* mtxmul(matrix *A, matrix *B)
{
    int Ar, Ac, Br, Bc;
    int i, j, k, l;
    matrix *C;

    C = NULL;
    
    Ar = nRows(A);
    Ac = nCols(A);
    Br = nRows(B);
    Bc = nCols(B);

    /* If the matricies dimensions aren't correct, return NULL */
    if(Ac != Br) {
        fprintf(stderr, "Error: Incompatible matrix dimensions.\n");
        return C;
    }
    
    /* Allocate Memory */
    C = CreateMatrix(Ar, Bc);
    
    /* Cik = AijBjk */
    for(i=0; i<Ar; i++) {
        for(k=0; k<Bc; k++) {
            for(j=0; j<Ac; j++) {
                C->array[i][k] += A->array[i][j] * B->array[j][k];
            }
        }
    }

    return C;
}

/* Code here courtesy of http://www.daniweb.com/software-development/c/code/216687 */
matrix* CalcMinor(matrix* A, int row, int col) {
    int i, j, a, b;
    int order;
    matrix *minor;

    minor = NULL;
    order = nRows(A);

    if(order <= 1)
        return NULL;
    if(row >= order || col >= order)
        return NULL;
    if( !(minor = CreateMatrix(order-1, order-1)) )
        return NULL;

    a = b = 0;

    for(i=0; i<order; i++) {
        if(i != row) {
            b = 0;
            for(j=0; j<order; j++) {
                if(j != col) {
                    setval(minor, val(A, i, j), a, b);
                    b++;
                }
            }
        a++;
        }
    }

    return minor;
}

/* Borrowed from the same site as CalcMinor */
double CalcDeterminant(matrix *p)
{
    int i, temp, order;
    double result;
    matrix *minor;

    result = 0;
    minor = NULL;
    order = nRows(p);

    if(order < 1) {
        fprintf(stderr, "CalcDeterminant(): Invalid Matrix.");
        return 0;
    }

    if(order == 1)
        return val(p, 0, 0);

    for(i=0; i<order; i++) {
        if( !(minor = CalcMinor(p, 0, i)) ) {
            fprintf(stderr, "CalcDeterminant(): Memory allocation failed.");
            return 0;
        }

        result += ( pow(-1, i) * p->array[0][i] * CalcDeterminant(minor));

        DestroyMatrix(minor);
    }

    return result;
}

matrix* CalcAdj(matrix* A)
{
    int i, j;
    double cofactor;
    matrix *minor, *adj, *adjt;
    
    adj = CreateMatrix(nRows(A), nRows(A));
    for(i=0; i<nRows(A); i++) {
        for(j=0; j<nRows(A); j++) {
            minor = CalcMinor(A, i, j);
            cofactor = pow(-1, (i+j+2)) * CalcDeterminant(minor);
            DestroyMatrix(minor);
            setval(adj, cofactor, i, j);
        }
    }

    adjt = mtxtrn(adj);
    DestroyMatrix(adj);
    

    return adjt;
}

matrix* CalcInv(matrix* A)
{
    matrix *inv, *adj;
    int i, j;
    double det;

    inv = NULL;
    det = CalcDeterminant(A);
    adj = CalcAdj(A);
    inv = CreateMatrix(nRows(A), nRows(A));
    
    for(i=0; i<nRows(A); i++) {
        for(j=0; j<nRows(A); j++) {
            setval(inv, val(adj, i, j)/det, i, j);
        }
    }

    return inv;
}

void Map(matrix* A, double (*func)(double))
{
    int i, j;
    for(i=0; i<nRows(A); i++) {
        for(j=0; j<nCols(A); j++) {
            setval(A, (*func)(val(A, i, j)), i, j);
        }
    }
}

/* Create a matrix of the specified dimensions and initialize all the values
 * and initialize all the values to zero.
 */
matrix* CreateMatrix(int row, int col)
{
    matrix *A;
    int i;
    A = NULL;

    if((row == 0) || (col == 0)) {
        printf("Matrix too small.");
        return A;
    }

    A = (matrix*) malloc(sizeof(matrix));
    if(!A) {
        fprintf(stderr, "Memory allocation error: %s\n", strerror(errno));
    }

    A->array = NULL;
    A->rows = 0;
    A->cols = 0;

    A->array = (double**) calloc(row, sizeof(double*));
    if(!A->array) {
        fprintf(stderr, "Memory allocation error: %s\n", strerror(errno));
        return A;
    }

    for(i=0; i<row; i++) {
        A->array[i] = (double*) calloc(col, sizeof(double));
        if(!A->array[i]) {
            fprintf(stderr, "Failed to allocate %lu bytes: %s\n", sizeof(double)*col, strerror(errno));
            return A;
        }
    }

    A->rows = row;
    A->cols = col;
    
    return A;
}

#define MAXROWS 100
#define MAXCOLS 100
#define LINELENGTH 80
matrix* ParseMatrix(char* raw)
{
    char *processed, *tmp;
    char **rows;
    double **values;
    matrix *out;
    int i, j;
    int nrows, ncols;

    processed = calloc(sizeof(char), MAXROWS*LINELENGTH);
    tmp = processed;

    while(*raw) {
        if(*raw == '[')
            *tmp = ' ';
        else if(*raw == ']')
            *tmp = ' ';
        else
            *tmp = *raw;
        tmp++;
        raw++;
    }

    /* Allocate memory to store the rows of the matrix as strings */
    values = (double**) calloc(sizeof(double*), MAXROWS);
    rows = (char**) malloc(sizeof(char*)*MAXROWS);
    for(i=0; i<MAXROWS; i++) {
        values[i] = (double*) calloc(sizeof(double), MAXCOLS);
        rows[i] = (char*) calloc(sizeof(char), LINELENGTH);
    }

    /* Row 1: */
    i = 0;
    tmp = strtok(processed, ";");
    strncpy(rows[i], tmp, LINELENGTH);

    /* All subsequent rows */
    for(i=1; tmp = strtok(NULL, ";"); i++) {
        strncpy(rows[i], tmp, LINELENGTH);
    }

    nrows = i;

    for(i=0; i<nrows; i++) {
        /* First column */
        j=0;
        tmp = strtok(rows[i], ",");
        values[i][j] = atof(tmp);

        ncols = 1;
        /* Rest of the columns */
        for(j=1; tmp = strtok(NULL, ","); j++) {
            if(j>ncols-1)
                ncols = j+1;
            values[i][j] = atof(tmp);
        }
    }

    out = CreateMatrix(nrows, ncols);
    for(i=0; i<nrows; i++) {
        for(j=0; j<ncols; j++)
            setval(out, values[i][j], i, j);
    }

    /* Free stuff */
    free(processed);
    for(i=0; i<MAXROWS; i++) {
        free(values[i]);
        free(rows[i]);
    }
    free(rows);
    free(values);
    
    return out;
}


/* Free the memory allocated by CreateMatrix */
void DestroyMatrix(matrix *A)
{
    int i;

    for(i=0; i<A->rows; i++) {
        free(A->array[i]);
    }
    free(A->array);
    free(A);
}

