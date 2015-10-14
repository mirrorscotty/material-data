#include "matrix.h"
#include "unifac.h"
#include <stdlib.h>
#include <stdio.h>

unifac_molec* UnifacCreateMolec(vector *ids, vector *count, unifac_data *d)
{
    int i, j;
    unifac_molec *m;

    m = (unifac_molec*) calloc(sizeof(unifac_molec), 1);
    m->ngroups = len(ids);
    m->ids = (int*) calloc(sizeof(int), m->ngroups);
    m->count = (int*) calloc(sizeof(int), m->ngroups);
    m->dat = (unifac_row*) calloc(sizeof(unifac_row), m->ngroups);

    for(i=0; i<m->ngroups; i++) {
        m->ids[i] = valV(ids, i);
        m->count[i] = valV(count, i);
        for(j=0; j<d->ngroups; j++) {
            if(d->rows[j].id == m->ids[i]) {
                m->dat[i] = d->rows[j];
                break;
            }
        }
    }

    return m;
}

void UnifacPrintMolecule(unifac_molec *m)
{
    int i;
    for(i=0; i<m->ngroups; i++)
        printf("\t%s\t(%d)\n", m->dat[i].subgroup, m->count[i]);
    return;
}
