#ifndef UNIFAC_H
#define UNIFAC_H

#include "matrix.h"

typedef struct {
    int id;
    int groupid;
    char group[20];
    char subgroup[20];
    double Rk;
    double Qk;
} unifac_row;

typedef struct {
    int ngroups;
    int* ids;
    int* count;
    unifac_row* dat;
} unifac_molec;

typedef struct {
    int ngroups;
    unifac_row *rows;
    matrix *interactions;
} unifac_data;

typedef struct {
    int nsolutes;
    double *xi;
    unifac_molec *m;
    unifac_data *dat;
} unifac_solution;

unifac_data* UnifacCreateData();
void UnifacPrintRow(unifac_row);
void UnifacPrintTable(unifac_data*);
void UnifacLoadData(unifac_data*, char*);

#endif

