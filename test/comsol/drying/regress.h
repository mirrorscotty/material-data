#ifndef REGRESS_H
#define REGRESS_H

/* Macro to shorten code where a number is added to the current value of an
 * element in a matrix. */
#define addval(A, VAL, I, J) setval((A), (VAL) + val((A), (I), (J)), (I), (J))

typedef struct {
    double **array;
    int rows;
    int cols;
} matrix;

void DestroyMatrix(matrix*);
matrix* CalcMinor(matrix*, int, int);
double CalcDeterminant(matrix*);
int nCols(matrix*);
int nRows(matrix*);
double val(matrix*, int, int);
void setval(matrix*, double, int, int);
matrix* CreateMatrix(int, int);
matrix* mtxtrn(matrix*);
matrix* mtxmul(matrix*, matrix*);
matrix* CalcAdj(matrix*);
matrix* CalcInv(matrix*);
matrix* regress(matrix*, matrix*);
matrix* polyfit(matrix*, matrix*, int);
matrix* ParseMatrix(char*);
matrix* linspace(double, double, int);
void Map(matrix*, double (*func)(double));
void mtxprnt(matrix*);
void mtxprntfile(matrix*, char*);

#endif
