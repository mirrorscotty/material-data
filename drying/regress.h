#ifndef REGRESS_H
#define REGRESS_H

typedef struct {
    double **array;
    int rows;
    int cols;
} matrix;

void DestroyMatrix(matrix*);
matrix* CalcMinor(matrix*, int, int);
double CalcDeterminant(matrix*);
int mtxlen1(matrix*);
int mtxlen2(matrix*);
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

#endif
