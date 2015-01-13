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
    double rhow, rhos, vs;
    choi_okos *co;

    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    vs = 1/(1+(rhos*Xdb/rhow));

    return (strain - vs)/strain;
}

