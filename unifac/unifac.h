typedef struct {
    int id;
    char[20] = group;
    char[20] = subgroup;
    double Rk;
    double Qk;
}, unifac_row;

typedef struct {
    int ngroups;
    int* ids;
    int* count;
    unifac_row* dat;
}, unifac_molec;

typedef struct {
    int ngroups;
    unifac_row *rows;
    matrix *interactions;
}, unifac_data;

typedef struct {
    int nsolutes;
    double *xi;
    unifac_molec *m;
    unifac_data *dat;
}, unifac_solution;

