#include <complex.h>
#include <math.h>
#include "matrix.h"

#define valVM(VECTOR, I) valV((VECTOR), (I)-1)
#define setvalVM(VECTOR, I, VAL) setvalV((VECTOR), (I)-1, (VAL))
#define valM(MATRIX, I, J) val((MATRIX), (I)-1, (J)-1)
#define setvalM(MATRIX, VAL, I, J) setval((MATRIX), (VAL), (I)-1, (J)-1)

double bnml(double n, double z)
{
    int i;
    double value = 1;
    for(i=1; i<=z; i++)
        value *= (n-z+i)/i;
    return value;
}

/**
 * Source: http://www.mathworks.com/matlabcentral/fileexchange/39035-numerical-inverse-laplace-transform/content//euler_inversion.m
 */
vector* ilt_euler(double complex (*f_s)(double complex), vector *t, int M)
{
    int i, j;
    double tmp, ilt;
    vector *xi, *k, *rebeta, *imbeta, *eta, *result;
    matrix *rebeta_mesh, *imbeta_mesh, *eta_mesh, *t_mesh;

    /* Set a default value for M if it isn't supplied */
    if(M < 1)
        M = 32;

    /* Make the xi vector */
    xi = CreateVector(2*M+1);
    setvalVM(xi, 1, .5);
    for(i=2; i<=M+1; i++)
        setvalVM(xi, i, 1);
    setvalV(xi, 2*M, pow(2, -1.*M));

    for(i=1; i<M; i++)
        setvalVM(xi,
                (2*M-i + 1),
                valVM(xi, 2*M-i+2) + pow(2, -1*M) * bnml(M, i));

    /* k is the iteration index */
    k = CreateVector(2*M+1);
    for(i=0; i<len(k); i++)
        setvalV(k, i, i);

    /* Create the real and imaginary parts of the beta vector separatedly. This
     * ensures that I don't need to make a complex vector and all the
     * operations for it. */
    rebeta = CreateVector(len(k));
    for(i=0; i<len(k); i++)
        setvalV(rebeta, i, M*log(10)/3);
    imbeta = CreateVector(len(k));
    for(i=0; i<len(k); i++)
        setvalV(imbeta, i, M_PI*valV(k, i));
    
    /* Create the eta vector */
    eta = CreateVector(len(k));
    for(i=0; i<len(eta); i++)
        setvalV(eta, i, (1-((int)valV(k, i) % 2)*2) * valV(xi, i));

    rebeta_mesh = meshgridX(rebeta, t);
    imbeta_mesh = meshgridX(imbeta, t);
    eta_mesh = meshgridX(eta, t);
    t_mesh = meshgridY(rebeta, t);

    /* Calculate the final result */
    result = CreateVector(len(t));
    for(i=0; i<len(t); i++) {
        ilt = pow(10, M/3.) / valV(t, i);
            
        tmp = 0;
        for(j=0; j<nCols(t_mesh); j++) {
            tmp += val(eta_mesh, i, j)
                *creal(f_s((val(rebeta_mesh, i, j)
                            + I*val(imbeta_mesh, i, j))/val(t_mesh, i, j)));
        }
        ilt *= tmp;
        setvalV(result, i, ilt);
    }

    DestroyVector(xi);
    DestroyVector(k);
    DestroyVector(rebeta);
    DestroyVector(imbeta);
    DestroyVector(eta);

    DestroyMatrix(rebeta_mesh);
    DestroyMatrix(imbeta_mesh);
    DestroyMatrix(eta_mesh);
    DestroyMatrix(t_mesh);

    return result;
}

