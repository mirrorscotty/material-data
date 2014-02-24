#ifndef DIFF_DATA_H
#define DIFF_DATA_H

struct _diff_data {
    /* Values */
    double kw; /* Permeability of water */
    double muw; /* Viscosity of water */
    double R; /* Ideal gas constant */
    double phi; /* Porosity */
};

/* Values for the Oswin isotherm correlation */
struct _oswin {
    double k0;
    double k1;
    double n0;
    double n1;
};
typedef struct _oswin oswin;

/* Values for the GAB isotherm equation */
struct _gab {
    double m0;
    double C0;
    double k0;
    double dHm;
    double dHk;
    double dHC;
};
typedef struct _gab gab;

oswin* CreateOswinData();
void DestroyOswinData(oswin*);
gab* CreateGABData();
void DestroyGABData();

#endif

