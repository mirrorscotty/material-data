#include "unifac.h"
#include <stdlib.h>
#include <stdio.h>

unifac_solution* UnifacCreateSolution(unifac_data *dat)
{
    unifac_solution *s;
    s = (unifac_solution*) calloc(sizeof(unifac_solution), 1);
    s->dat = dat;
    return s;
}

void UnifacDestroySolution(unifac_solution *s)
{
    free(s);
}

void UnifacAddMolec(unifac_solution *s, unifac_molec *m, double molefrac)
{
    double *xinew;
    int i;
    unifac_molec *molecnew;

    xinew = (double*) calloc(sizeof(double), s->nsolutes + 1);
    molecnew = (unifac_molec*) calloc(sizeof(unifac_molec), s->nsolutes + 1);

    if(s->xi) {
        for(i=0; i<s->nsolutes; i++) {
            xinew[i] = s->xi[i];
            molecnew[i] = s->m[i];
        }
        xinew[s->nsolutes] = molefrac;
        molecnew[s->nsolutes] = *m;
    } else {
        xinew[0] = molefrac;
        molecnew[0] = *m;
    }
    s->nsolutes = s->nsolutes + 1;
    free(s->xi);
    free(s->m);
    s->xi = xinew;
    s->m = molecnew;
    return;
}

unifac_solution* UnifacPureSolution(int i, unifac_solution *s)
{
    unifac_solution *pure;
    unifac_molec *m;
    double xi = 1;

    m = &(s->m[i]);

    pure = UnifacCreateSolution(s->dat);
    UnifacAddMolec(pure, m, xi);

    return pure;
}

void UnifacPrintSolution(unifac_solution *s)
{
    int i;
    printf("%d Solutes\n", s->nsolutes);
    for(i=0; i<s->nsolutes; i++) {
        printf("Molecule #%d (x = %g):\n", i, s->xi[i]);
        UnifacPrintMolecule( &(s->m[i]) );
    }
}

