#include "pasta.h"
#include "choi-okos.h"

/**
 * Effective density
 * @param cw Mass concentration of water [kg/m^3]
 * @param wv Mass fraction of water vapor [-]
 * @param phi porosity [-]
 * @param T temperature [K]
 * @param P pressure [Pa]
 * @returns Density [kg/m^3]
 */
double rho_eff(double cw, double wv, double phi, double T, double P)
{
    double rho_eff, /* Effective density */
           rhos, /* Solid density */
           rhow, /* Liquid density (water) */
           rhog, /* Gas density */
           Sw, /* Water saturation */
           Sg; /* Gas saturation */

    choi_okos *co;
    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    rhog = rho_gas(wv, T, P);
    Sw = sat_wat(cw, phi, T);
    Sg = sat_gas(cw, phi, T);
    
    rho_eff = (1-phi)*rhos + phi*(Sw*rhow + Sg*rhog);
    return rho_eff;
}

/**
 * Effective heat capacity
 * @param cw mass concentration of water [kg/m^3]
 * @param wv mass fraction of water vapor [-]
 * @param phi porosity [-]
 * @param T temperature [K]
 * @param P pressure [Pa]
 * @returns Specific heat capacity [J/K]
 */
double Cp_eff(double cw, double wv, double phi, double T, double P)
{
    double Cp_eff,
           Mg, /* Mass fraction gas */
           Mw, /* Mass fractiong water */
           Ms, /* Mass fraction solid */
           wa = 1-wv, /* Mass fraction of air in gas phase */
           Cpv = CPVAPOR, /* vapor heat capacity  at 100 C [J/(g K)] */
           Cpa = CPAIR, /* air heat capacity at 23 C, 760 mmHg [J/(g K) */
           Cpw, /* Water heat capacity */
           Cps, /* Solid heat capacity */
           rhos, /* solid density */
           cs,
           cg;
    choi_okos *co;
    co = CreateChoiOkos(WATERCOMP);
    Cpw = Cp(co, T);
    DestroyChoiOkos(co);

    co = CreateChoiOkos(PASTACOMP);
    Cps = Cp(co, T);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    /* Calculate mass fractions */
    cs = rhos*phi; /* Solid mass concentration */
    cg = conc_gas(cw, wv, phi, T, P); /* Gas phase mass concentration */
    Mg = cg/(cg+cw+cs); /* Mass fraction gas */
    Mw = cw/(cg+cw+cs); /* Mass fraction water */
    Ms = cs/(cg+cw+cs); /* Mass fraction solid */
    
    Cp_eff = Mg*(wv*Cpv + wa*Cpa) + Mw*Cpw + Ms*Cps;
    return Cp_eff;
}

/**
 * Effective thermal conductivity
 * @param cw mass concentration of water [kg/m^3]
 * @param wv mass fraction of water vapor [-]
 * @param phi porosity [-]
 * @param T temperature [K]
 * @param P pressure [Pa]
 * @returns Thermal conductivity [W/m-K]
 */
double k_eff(double cw, double wv, double phi, double T, double P)
{
    double k_eff,
           ks, /* solid thermal conductivity */
           kw, /* k for liquid water */
           ka = KAIR,
           kv = KVAPOR,
           wa = 1-wv, /* Mass fraction of air in gas phase */
           Sw, /* Liquid water saturation */
           Sg; /* Gas saturation */
    choi_okos *co;

    co = CreateChoiOkos(WATERCOMP);
    kw = k(co, T);
    DestroyChoiOkos(co);

    co = CreateChoiOkos(PASTACOMP);
    ks = k(co, T);
    DestroyChoiOkos(co);

    Sw = sat_wat(cw, phi, T);
    Sg = sat_gas(cw, phi, T);

    k_eff = (1-phi)*ks + phi*(Sw*kw + Sg*(wa*ka + wv*kv));
    return k_eff;
}

/* Heat flux due to fluid (liquid and gas) flow
 * TODO Add in velocities
 * cw: mass concentration of water [kg/m^3]
 * wv: mass fraction of water vapor
 * phi: porosity
 * T: temperature [K]
 * P: pressure [Pa]
 */
double fluidconv(double cw, double wv, double phi, double T, double P)
{
    double rhoCpv,
           rhow, /* density of water phase */
           rhog = rho_gas(wv, T, P), /* density of gas phase */
           vw = 0, /* water phase velocity */
           vg = 0, /* gas phase velocity */
           Cpw, /* water heat capacity */
           Cpv = CPVAPOR,
           Cpa = CPAIR,
           wa = 1-wv, /* air mass fraction in gas */
           gradcw = 0, /* Concentration gradient */
           D = 0; /* Capillary diffusivity */
    choi_okos *co;

    co = CreateChoiOkos(WATERCOMP);
    Cpw = Cp(co, T);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    rhoCpv = (rhow*vw - D*gradcw)*Cpw + rhog*vg*(wv*Cpv+wa*Cpa);
    return rhoCpv;
}

