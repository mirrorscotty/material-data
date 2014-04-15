#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include "pasta.h"
#include "isotherms.h"
#include "choi-okos.h"

#define BUFFSIZE 80
#define OUTSTR "%s,%g,%g,%g,%g\n"

char* fivevarcw(const char *fname,
        double (*f)(double, double, double, double, double),
        double del,
        double cw,
        double wv,
        double phi,
        double T,
        double P)
{
    double low, mid, high, change;
    char *str;

    str = (char*) calloc(sizeof(char), BUFFSIZE);
    /* Name, Low, Mid, High, %Change */
    low = f(cw*(1-del), wv, phi, T, P);
    mid = f(cw, wv, phi, T, P);
    high = f(cw*(1+del), wv, phi, T, P);
    change = (high-low)/low;

    sprintf(str, OUTSTR, fname, low, mid, high, change);

    return str;
}

char* fivevarwv(const char *fname,
        double (*f)(double, double, double, double, double),
        double del,
        double cw,
        double wv,
        double phi,
        double T,
        double P)
{
    double low, mid, high, change;
    char *str;

    str = (char*) calloc(sizeof(char), BUFFSIZE);
    /* Name, Low, Mid, High, %Change */
    low = f(cw, wv*(1-del), phi, T, P);
    mid = f(cw, wv, phi, T, P);
    high = f(cw, wv*(1+del), phi, T, P);
    change = (high-low)/low;

    sprintf(str, OUTSTR, fname, low, mid, high, change);

    return str;
}

char* fivevarphi(const char *fname,
        double (*f)(double, double, double, double, double),
        double del,
        double cw,
        double wv,
        double phi,
        double T,
        double P)
{
    double low, mid, high, change;
    char *str;

    str = (char*) calloc(sizeof(char), BUFFSIZE);
    /* Name, Low, Mid, High, %Change */
    low = f(cw, wv, phi*(1-del), T, P);
    mid = f(cw, wv, phi, T, P);
    high = f(cw, wv, phi*(1+del), T, P);
    change = (high-low)/low;

    sprintf(str, OUTSTR, fname, low, mid, high, change);

    return str;
}

char* fivevarT(const char *fname,
        double (*f)(double, double, double, double, double),
        double del,
        double cw,
        double wv,
        double phi,
        double T,
        double P)
{
    double low, mid, high, change;
    char *str;

    str = (char*) calloc(sizeof(char), BUFFSIZE);
    /* Name, Low, Mid, High, %Change */
    low = f(cw, wv, phi, T*(1-del), P);
    mid = f(cw, wv, phi, T, P);
    high = f(cw, wv, phi, T*(1+del), P);
    change = (high-low)/low;

    sprintf(str, OUTSTR, fname, low, mid, high, change);

    return str;
}

char* fivevarP(const char *fname,
        double (*f)(double, double, double, double, double),
        double del,
        double cw,
        double wv,
        double phi,
        double T,
        double P)
{
    double low, mid, high, change;
    char *str;

    str = (char*) calloc(sizeof(char), BUFFSIZE);
    /* Name, Low, Mid, High, %Change */
    low = f(cw, wv, phi, T, P*(1-del));
    mid = f(cw, wv, phi, T, P);
    high = f(cw, wv, phi, T, P*(1+del));
    change = (high-low)/low;

    sprintf(str, OUTSTR, fname, low, mid, high, change);

    return str;
}

char* threevarcw(const char *fname,
        double (*f)(double, double, double),
        double del,
        double cw,
        double phi,
        double T)
{
    double low, mid, high, change;
    char *str;

    str = (char*) calloc(sizeof(char), BUFFSIZE);
    /* Name, Low, Mid, High, %Change */
    low = f(cw*(1-del), phi, T);
    mid = f(cw, phi, T);
    high = f(cw*(1+del), phi, T);
    change = (high-low)/low;

    sprintf(str, OUTSTR, fname, low, mid, high, change);

    return str;
}

char* threevarphi(const char *fname,
        double (*f)(double, double, double),
        double del,
        double cw,
        double phi,
        double T)
{
    double low, mid, high, change;
    char *str;

    str = (char*) calloc(sizeof(char), BUFFSIZE);
    /* Name, Low, Mid, High, %Change */
    low = f(cw, phi*(1-del), T);
    mid = f(cw, phi, T);
    high = f(cw, phi*(1+del), T);
    change = (high-low)/low;

    sprintf(str, OUTSTR, fname, low, mid, high, change);

    return str;
}

char* threevarT(const char *fname,
        double (*f)(double, double, double),
        double del,
        double cw,
        double phi,
        double T)
{
    double low, mid, high, change;
    char *str;

    str = (char*) calloc(sizeof(char), BUFFSIZE);
    /* Name, Low, Mid, High, %Change */
    low = f(cw, phi, T*(1-del));
    mid = f(cw, phi, T);
    high = f(cw, phi, T*(1+del));
    change = (high-low)/low;

    sprintf(str, OUTSTR, fname, low, mid, high, change);

    return str;
}

char* onevar(const char *fname,
        double (*f)(double),
        double del,
        double T)
{
    double low, mid, high, change;
    char *str;

    str = (char*) calloc(sizeof(char), BUFFSIZE);
    /* Name, Low, Mid, High, %Change */
    low = f(T*(1-del));
    mid = f(T);
    high = f(T*(1+del));
    change = (high-low)/low;

    sprintf(str, OUTSTR, fname, low, mid, high, change);

    return str;
}

/* Print out stuffs */
int main(int argc, char *argv[])
{
    /* Make a new directory to store all the stuff in if it doesn't exist. */
    struct stat st = {0};
    if (stat("sensitivity", &st) == -1)
        mkdir("sensitivity", 0700);

    double cw,
           wv = 0.5-.25,
           phi = .09/2,
           T = (60/2)+273.15,
           P = 10130/2,
           del = .1;
    double Xdb,
           rhos;
    choi_okos *co;
    oswin *o;
    FILE *fp;

    Xdb = .1 * mdb_wat_sat(phi, T);

    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    cw = Xdb*rhos*(1-phi);

    /* cw */
    fp = fopen("sensitivity/cw.csv", "w");
    fprintf(fp, "Sensitivity Analysis\nwv = %g\nphi = %g\nT = %g\nP = %g\n",
            wv, phi, T, P);
    fprintf(fp, "cw,%g,%g,%g,Change\n",cw*(1-del), cw, cw*(1+del));
    fprintf(fp, "%s", fivevarcw("cg", &conc_gas, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarcw("ca", &conc_air, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarcw("cv", &conc_vap, del, cw, wv, phi, T, P));

    fprintf(fp, "%s", fivevarcw("rho_eff", &rho_eff, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarcw("Cp_eff", &Cp_eff, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarcw("k_eff", &k_eff, del, cw, wv, phi, T, P));

    fprintf(fp, "%s", threevarcw("Xdb", &mdb_wat, del, cw, phi, T));
    fprintf(fp, "%s", threevarcw("Sw", &sat_wat, del, cw, phi, T));
    fprintf(fp, "%s", threevarcw("Sg", &sat_gas, del, cw, phi, T));
    fprintf(fp, "%s", threevarcw("kw", &perm_wat, del, cw, phi, T));
    fclose(fp);

    /* wv */
    fp = fopen("sensitivity/wv.csv", "w");
    fprintf(fp, "Sensitivity Analysis\ncw = %g\nphi = %g\nT = %g\nP = %g\n",
            cw, phi, T, P);
    fprintf(fp, "wv,%g,%g,%g,Change\n",wv*(1-del), wv, wv*(1+del));
    fprintf(fp, "%s", fivevarwv("cg", &conc_gas, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarwv("ca", &conc_air, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarwv("cv", &conc_vap, del, cw, wv, phi, T, P));

    fprintf(fp, "%s", fivevarwv("rho_eff", &rho_eff, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarwv("Cp_eff", &Cp_eff, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarwv("k_eff", &k_eff, del, cw, wv, phi, T, P));

    fprintf(fp, "%s", threevarcw("rho_g", &rho_gas, del, wv, T, P));

    fprintf(fp, "%s", onevar("xv", &molefrac_vap, del, wv));
    fclose(fp);

    /* phi */
    fp = fopen("sensitivity/phi.csv", "w");
    fprintf(fp, "Sensitivity Analysis\ncw = %g\nwv = %g\nT = %g\nP = %g\n",
            cw, phi, T, P);
    fprintf(fp, "phi,%g,%g,%g,Change\n",phi*(1-del), phi, phi*(1+del));
    fprintf(fp, "%s", fivevarphi("cg", &conc_gas, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarphi("ca", &conc_air, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarphi("cv", &conc_vap, del, cw, wv, phi, T, P));

    fprintf(fp, "%s", fivevarphi("rho_eff", &rho_eff, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarphi("Cp_eff", &Cp_eff, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarphi("k_eff", &k_eff, del, cw, wv, phi, T, P));

    fprintf(fp, "%s", threevarphi("Xdb", &mdb_wat, del, cw, phi, T));
    fprintf(fp, "%s", threevarphi("Sw", &sat_wat, del, cw, phi, T));
    fprintf(fp, "%s", threevarphi("Sg", &sat_gas, del, cw, phi, T));
    fprintf(fp, "%s", threevarphi("kw", &perm_wat, del, cw, phi, T));
    fclose(fp);

    /* T */
    fp = fopen("sensitivity/T.csv", "w");
    fprintf(fp, "Sensitivity Analysis\ncw = %g\nwv = %g\nphi = %g\nP = %g\n",
            cw, wv, phi, P);
    fprintf(fp, "T,%g,%g,%g,Change\n",T*(1-del), T, T*(1+del));
    fprintf(fp, "%s", fivevarT("cg", &conc_gas, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarT("ca", &conc_air, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarT("cv", &conc_vap, del, cw, wv, phi, T, P));

    fprintf(fp, "%s", fivevarT("rho_eff", &rho_eff, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarT("Cp_eff", &Cp_eff, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarT("k_eff", &k_eff, del, cw, wv, phi, T, P));

    fprintf(fp, "%s", threevarT("Xdb", &mdb_wat, del, cw, phi, T));
    fprintf(fp, "%s", threevarT("Sw", &sat_wat, del, cw, phi, T));
    fprintf(fp, "%s", threevarT("Sg", &sat_gas, del, cw, phi, T));
    fprintf(fp, "%s", threevarT("kw", &perm_wat, del, cw, phi, T));

    fprintf(fp, "%s", threevarT("rho_g", &rho_gas, del, wv, phi, T));

    fprintf(fp, "%s", onevar("mu_w", &visc_wat, del, T));
    fprintf(fp, "%s", onevar("Pv", &pvap_wat, del, T));
    fclose(fp);

    /* P */
    fp = fopen("sensitivity/P.csv", "w");
    fprintf(fp, "Sensitivity Analysis\ncw = %g\nwv = %g\nphi = %g\nT = %g\n",
            cw, phi, T, P);
    fprintf(fp, "P,%g,%g,%g,Change\n",P*(1-del), P, P*(1+del));
    fprintf(fp, "%s", fivevarP("cg", &conc_gas, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarP("ca", &conc_air, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarP("cv", &conc_vap, del, cw, wv, phi, T, P));

    fprintf(fp, "%s", fivevarP("rho_eff", &rho_eff, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarP("Cp_eff", &Cp_eff, del, cw, wv, phi, T, P));
    fprintf(fp, "%s", fivevarP("k_eff", &k_eff, del, cw, wv, phi, T, P));

    fprintf(fp, "%s", threevarT("rho_g", &rho_gas, del, wv, T, P));
    fclose(fp);

    return 0;
}

