#include "material-data.h"

/**
 * Calculate porosity based on moisture content and volumetric strain. Here,
 * porosity is defined as the portion of the volume not taken up by the solid
 * phase. This definition is useful for calculating such things as effective
 * stress and Poisson ratio. The equation here was based on the assumptions
 * that:
 * 1. Adsorbed water has the same density as bulk water.
 * 2. The pasta sample initially contains no air.
 * 3. The volume of the solid phase is constant.
 *
 * \f[
 * \phi = \frac{\epsilon - V_s}{\epsilon}
 * \f]
 * \f[
 * V_s = \frac{1}{1+\frac{\rho_s X_{db}}{\rho_w}}
 * \f]
 *
 * @param Xdb Dry basis moisture content [kg/kg db]
 * @param T Temperature [K]
 * @param strain Volumetric strain [-]
 * @returns Porosity [-]
 */
double porosity(double Xdb, double T, double strain)
{
    double rhow, rhos, vs, vw, vo, phi;
    choi_okos *co;

    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    vs = 1/rhos;
    vw = Xdb/rhow;
    vo = vs + vw;

    phi = strain - (vs-vo)/vo;
    if(phi < 0)
        phi = 0;
    return phi;
}

/** Unlike the porosity() function, this should calculate solid fraction
 * correctly. */
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

