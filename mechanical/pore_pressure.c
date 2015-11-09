#include "mechanical.h"
#include "pasta.h"
#include "choi-okos.h"
#include "glass-transition.h"
#include <math.h>
#define SSP .1

/**
 * Calculate pore pressure based on the Kelvin equation.
 * @param Xdb Moisture content [kg/kg db]
 * @param T Temperature [K]
 * @returns Capillary pressure [Pa]
 */
double pore_press(double Xdb, double T)
{
    double rhow,
           R = GASCONST,
           Vm = 1.802e-5, /* m^3/mol */
           aw;
    choi_okos *co;
    oswin *o;
//    gordontaylor *gt;
//    gt = GTSemolina();

    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    o = CreateOswinData();
    aw = OswinInverse(o, Xdb, T);
    DestroyOswinData(o);

    /* If there isn't enough water to form a meniscus, then there is no
     * capillary pressure. The cutoff here is completely made up. */
    //if(Xdb > GordonTaylorInv(gt,T))
    //if(aw > .7)
        return R*T/Vm * log(aw);
    //else
    //    return 0;
}

