#include <stdio.h>
#include <stdlib.h>
#include "regress.h"

/* Equation: u'' = Pe*u' */
/* Create the element matrix for the specified values of Pe and h */
matrix* CreateElementMatrix(double Pe, double h)
{
    matrix *elem;
 
    elem = CreateMatrix(2, 2);
    setval(elem, 1/h-Pe/2, 0, 0);
    setval(elem, -1/h+Pe/2, 0, 1);
    setval(elem, -1/h-Pe/2, 1, 0);
    setval(elem, 1/h+Pe/2, 1, 1);

    return elem;
}

/* Create the load vector */
matrix* CreateElementLoad(double Pe, double h) {
    matrix *f;

    f = CreateMatrix(2, 1);

    setval(f, 0, 0, 0);
    setval(f, 0, 1, 0);

    return f;
}

/* Equation: u'' + u = 1 */
/* Test functions to create the element and load matricies based on the equation
 * in Problem 1
 */
matrix* testelem(double Pe, double h)
{
    matrix *elem;
 
    elem = CreateMatrix(2, 2);
    setval(elem, -1/h+h/3, 0, 0);
    setval(elem, 1/h+h/6, 0, 1);
    setval(elem, 1/h+h/6, 1, 0);
    setval(elem, -1/h+h/3, 1, 1);

    return elem;
}

matrix* testload(double Pe, double h) {
    matrix *f;

    f = CreateMatrix(2, 1);

    setval(f, h/2, 0, 0);
    setval(f, h/2, 1, 0);

    return f;
}

/* Assemble the global coefficient matrix for a non-uniform mesh. The first
 * argument is a function pointer to the function which creates the element
 * matrix, and the second is the mesh to be used.
 */
matrix* AssembleJ(matrix* (*makej)(double, double), matrix* mesh, double Pe)
{
    matrix *J, *j;
    int n, i;

    /* Determine the number of elements from the mesh */
    n = mtxlen2(mesh);

    /* Create a blank global matrix */
    J = CreateMatrix(n+1, n+1);

    for(i=0; i<n; i++) {
	/* Generate the element matrix for the specified element width */
	j = makej(Pe, val(mesh, i, 0));

	/* Add the values of the element matrix to the global matrix */
	addval(J, val(j, 0, 0), i, i);
	addval(J, val(j, 0, 1), i, i+1);
	addval(J, val(j, 1, 0), i+1, i);
	addval(J, val(j, 1, 1), i+1, i+1);

	/* Clean up */
	DestroyMatrix(j);
    }

    return J;
}

matrix* AssembleF(matrix* (*makef)(double, double), matrix* mesh, double Pe)
{
    matrix *F, *f;
    int n, i;

    n = mtxlen2(mesh);

    F = CreateMatrix(n+1, 1);

    for(i=0; i<n; i++) {
	f = makef(Pe, val(mesh, i, 0));

	addval(F, val(f, 0, 0), i, 0);
	addval(F, val(f, 1, 0), i+1, 0);

	DestroyMatrix(f);
    }

    return F;
}


/* Simple function to alter J and F so that Dirichlet boundary conditions are
 * imposed on both ends of the domain. In this case, "leftbc" is imposed at
 * c = 0, and rightbc is imposed at c = 1.
 */
void ApplyBoundaryConditions(matrix* J, matrix* F, double leftbc, double rightbc)
{
    int rows = mtxlen2(J);
    int i;

    for(i=0;i<rows; i++) {
	setval(J, 0, 0, i);
	setval(J, 0, rows-1, i);
    }

    setval(J, 1, 0, 0);
    setval(J, 1, rows-1, rows-1);

    setval(F, leftbc, 0, 0);
    setval(F, rightbc, rows-1, 0);
}

/* Determine the element width based on the size of the domain and the number
 * of elements.
 */
matrix* GenerateUniformMesh(double left, double right, int n)
{
    double h;
    int i;
    matrix *mesh;
    mesh = CreateMatrix(n, 1);

    h = (right-left)/n;

    for(i=0; i<n; i++)
	setval(mesh, h, i, 0);

    return mesh;
}

matrix* MeshXCoords(matrix *mesh, double left, double right)
{
    matrix *x;
    int i;
    int n = mtxlen2(mesh);

    x = CreateMatrix(n+1, 1);

    setval(x, left, 0, 0);

    for(i=1; i<=n; i++) {
        setval(x, val(x, i-1, 0) + val(mesh, i-1, 0), i, 0);
    }

    if((val(x, n, 0) - right) > 1e-5) {
        fprintf(stderr, "Warning: mesh not aligned.\n");
    }

    return x;
}

#define PI 3.141592654

/* This program takes up to three arguments. The first is the Peclet number to
 * use in the calculation. The second is the number of elements to use when
 * solving the problem, and the third (optional) argument is the filename to
 * save the solution to. */
int main(int argc, char *argv[])
{
    double Pe, h;
    double left = 0;
    double right = 1;
    //double right = PI/2;
    int n;
    matrix *J, *F, *u, *Jinv, *mesh, *x;

    /* Parse arguments */
    if(argc < 2) {
	fprintf(stderr, "Too few arguments: exiting.\n");
	return 1;
    }
    Pe = atof(argv[1]);

    n = atoi(argv[2]);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh(left, right, n);

    /* Use a mesh with these node spacings */
    //mesh = ParseMatrix("[.75;.15;.0125;.0125;.0125;.0125;.0125;.0125;.0125;.0125]");

    x = MeshXCoords(mesh, left, right);

    /* All the nodes are equally spaced. */
    h = val(mesh, 0, 0);

    J = AssembleJ(&CreateElementMatrix, mesh, Pe);
    //J = AssembleJMesh(&testelem, mesh, Pe);

    F = AssembleF(&CreateElementLoad, mesh, Pe);
    //F = AssembleFMesh(&testload, mesh, Pe);

    ApplyBoundaryConditions(J, F, 0, 1);
    //ApplyBoundaryConditions(J, F, 2, 1);

    /* Invert the coefficient matrix and premultiply the load vector by it. */
    Jinv = CalcInv(J);
    u = mtxmul(Jinv, F);

    /* Print out the result */
    if(argc == 4) {
	mtxprntfile(u, argv[3]);
    } else {
        mtxprnt(x);
        puts("");
	mtxprnt(u);
    }
    
    /* Clean up the allocated memory */
    DestroyMatrix(J);
    DestroyMatrix(F);
    DestroyMatrix(u);
    DestroyMatrix(Jinv);
    DestroyMatrix(mesh);

    return 0;
}
