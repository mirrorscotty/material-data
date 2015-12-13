#include "material-data.h"

double solidfrac(double Xo, double T, double strain)
{
    double rhow, rhos, vv0, xf;
    choi_okos *co;

    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    xf = Xo * (rhos/rhow)*(1-strain) - strain;
    if(xf < 0)
        return 0;
    else if(xf > 1)
        return 1;
    else
        return xf;
} 

/**
 * Calculate porosity based on moisture content and volumetric strain. Here,
 * porosity is defined as the portion of the volume taken up by the air 
 * phase. The equation here was based on the assumptions
 * that:
 * 1. Adsorbed water has the same density as bulk water.
 * 2. The pasta sample initially contains no air.
 * 3. The volume of the solid phase is constant.
 *
 * @param Xo Initial moisture content [kg/kg db]
 * @param Xdb Dry basis moisture content [kg/kg db]
 * @param T Temperature [K]
 * @param strain Volumetric strain [-]
 * @returns Porosity [-]
 */
double porosity(double Xo, double Xdb, double T, double strain)
{
    double rhow, rhos, vv0, phi;
    choi_okos *co;

    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    phi = 1 - 1/(strain+1) * (rhos/rhow*Xdb + 1)/(rhos/rhow*Xo + 1);
    if(phi < 0)
        return 0;
    else if(phi > 1)
        return 1;
    else
        return phi;
} 

double DporosityDX(double Xo, double Xdb, double T, double strain)
{
    double rhow, rhos, vv0, phi, Dphi, rhoR;
    choi_okos *co;

    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    rhoR = rhos/rhow;

    phi = 1 - 1/(strain+1) * (rhos/rhow*Xdb + 1)/(rhos/rhow*Xo + 1);
    Dphi = -1*rhoR/(Xo*strain*rhoR + Xo*rhoR + strain + 1);
    
    if(phi < 0)
        return 0;
    else if(phi > 1)
        return 0;
    else
        return Dphi;
} 

double DporosityDstrain(double Xo, double Xdb, double T, double strain)
{
    double rhow, rhos, vv0, phi, Dphi, rhoR;
    choi_okos *co;

    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    rhoR = rhos/rhow;

    phi = 1 - 1/(strain+1) * (rhos/rhow*Xdb + 1)/(rhos/rhow*Xo + 1);
    Dphi = (Xdb*rhoR+1)/pow(strain+1,2)*(Xo*rhoR+1);
    
    if(phi < 0)
        return 0;
    else if(phi > 1)
        return 0;
    else
        return Dphi;
} 

