/**
 * @file unifac.c
 * Implementation of the UNIFAC method for predicting VLE.
 */

#include <math.h>
#include "unifac.h"

/* Combinatorial */

#define Z 10 /* Coordination number */

/**
 * The r_i value calculated from group surface area.
 * \f[
 * r_i = \sum_{k=1}^n v_k R_k
 * \f]
 * @param molec Molecule to calculate the surface area for.
 * @returns r_i
 */
double _r(unifac_molec molec)
{
    int k, sum = 0;
    for(k=0; k<molec.ngroups; k++)
        sum += molec.ngroups * molec.dat->Rk;
    return sum;
}

/**
 * The q_i value calculated from group volume.
 * \f[
 * q_i = \sum_{k=1}^n v_k Q_k
 * \f]
 * @param molec Molecule to calculate the volume for.
 * @returns q_i
 */
double _q(unifac_molec molec)
{
    int k, sum = 0;
    for(k=0; k<molec.ngroups; k++)
        sum += molec.ngroups * molec.dat->Qk;
    return sum;
}

/**
 * The L_i value is calculated from q, r, and z, where q and r are defined
 * above and z is the coordination number. Here, z is assumed to be equal to
 * ten.
 * \f[
 * L_i = \frac{z}{2}(r_i - q_i) - (r_i -1)
 * \f]
 * @param molec Molecule to calculate the value for
 * @returns L_i
 */
double _L(unifac_molec molec)
{
    double ri, qi;
    ri = _r(molec);
    qi = _q(molec);
    return Z/2*(ri - qi) - (ri - 1);
}

/**
 * Area fractional component for the ith molecule in the solution.
 * \f[
 * \phi_i = \frac{x_i r_i}{\sum_{j=1}^n x_j r_j}
 * \f]
 * @param i Molecule ID in the solution.
 * @param s Solution
 * @returns phi_i
 */
double _phi(int i, unifac_solution *s)
{
    int j;
    double numer, denom = 0;
    numer = s->xi[i] * _r(s->m[i]);
    for(j=0; j<s->nsolutes; j++)
        denom += s->xi[j] * _r(s->m[j]);
    return numer/denom;
}

double _phiElec(int i, unifac_solution *s)
{
    int j;
    double numer, denom = 0;
    numer = s->xi[i] * pow(_r(s->m[i]), 3./4.);
    for(j=0; j<s->nsolutes; j++)
        denom += s->xi[j] * pow(_r(s->m[j]), 3./4.);
    return numer/denom;
}

/**
 * Molar weighted segment for the ith molecule in the solution.
 * \f[
 * \theta_i = \frac{x_i q_i}{\sum_{j=1}^n x_j q_j}
 * \f]
 * @param i Molecule ID in the solution.
 * @param s Solution
 * @returns theta_i
 */
double _theta(int i, unifac_solution *s)
{
    int j;
    double numer, denom = 0;
    numer = s->xi[i] * _q(s->m[i]);
    for(j=0; j<s->nsolutes; j++)
        denom += s->xi[j] * _q(s->m[j]);
    return numer/denom;
}

/**
 * Calculate that natural log of the combinatorial component of the activity coefficient for the ith molecule in the solution.
 * \f[
 * \ln\gamma_i^c = \ln\frac{\phi_i}{x_i}
 *      + \frac{z}{2}q_i\ln\frac{\theta_i}{\phi_i} + L_i
 *      - \frac{\phi_i}{x_i}\sum_{j=1}^n x_j L_j
 * \f]
 * @param i Molecule ID in the solution.
 * @param s Solution
 * @returns ln(gamma_i^c)
 */
double _ln_gammac(int i, unifac_solution *s)
{
    int j;
    double result, sum = 0;

    for(j=0; j<s->nsolutes; j++)
        sum += s->xi[i] * _L(s->m[i]);

    result = log(_phi(i, s)/s->xi[i])
            + Z/2*_q(s->m[i])*log(_theta(i, s)/_phi(i, s))
            + _L(s->m[i])
            - _phi(i, s)/s->xi[i] * sum;

    return result;
}

double _ln_gammacElec(int i, unifac_solution *s)
{
    double result;

    result = log(_phiElec(i, s)/s->xi[i]) + 1 - _phiElec(i, s)/s->xi[i]
            + Z/2*_q(s->m[i]) * (1 - _phi(i, s)/_theta(i, s)
                                 + log(_phi(i, s)/_theta(i, s)));

    return result;
}

/**
 * Take a count of the number of groups with a particular ID in a specific
 * molecule. Return zero if the group isn't in there.
 * @param id Group ID
 * @param m Molecule
 * @returns Number of groups
 */
int molec_group_count(int id, unifac_molec m)
{
    int i = 0;
    while(i<m.ngroups) {
        if(id == m.ids[i])
            return m.count[i];
        i++;
    }
    return 0;
}

int solution_group_count(int id, unifac_solution *s)
{
    int i, sum = 0;
    for(i=0; i<s->nsolutes; i++)
        sum += molec_group_count(id, s->m[i]);
    return sum;
}


/* Residual */
/**
 * Group mole fraction. (Number of moles of the selected group divided by
 * total number of moles of groups.)
 * \f[
 * X_m = \frac{\sum_j v_m^j x_j}{\sum_j\sum_n v_n^j x_j}
 * \f]
 * @param m Group ID
 * @param s Solution
 * @returns X_m
 */
double _Xm(int m, unifac_solution *s)
{
    int j, n;
    double numer = 0, denom = 0, result;
    for(j=0; j<s->nsolutes; j++) {
        numer += molec_group_count(m, s->m[j]) * s->xi[j];
        for(n=0; n<s->dat->ngroups; n++)
            denom += molec_group_count(n, s->m[j]) * s->xi[j];
    }
    
    result = numer/denom;
    return result;
}

/**
 * Group interaction parameter. This quantifies the interaction energy
 * between groups. The units are in terms of Kelvins.
 * @param m First group ID
 * @param n Second group ID
 * @param s Solution
 * @param T Temperature [K]
 * @return Psi_mn
 */
double _Psi(int m, int n, unifac_solution *s, double T)
{
    int id1 = s->m[m].dat->id,
        id2 = s->m[n].dat->id,
        i;
    matrix *a;
    double a_mn;

    a = s->dat->interactions;

    for(i=0; i<nRows(a); i++) {
        if( (val(a, i, 0) == id1) && (val(a, i, 1) == id2) ) {
            a_mn = val(a, i, 2);
        }
    }

    return -a_mn/T;
}

double _PsiElec(int m, int n, unifac_solution *s, double T)
{
    int i;
    matrix *a;
    double a_mn, b_mn, c_mn;
    double result;

    if(solution_group_count(m, s) == 0)
        return 0;
    if(solution_group_count(n, s) == 0)
        return 0;

    a = s->dat->interactions;

    for(i=0; i<nRows(a); i++) {
        if( (val(a, i, 0) == m) && (val(a, i, 1) == n) ) {
            a_mn = val(a, i, 2);
            b_mn = val(a, i, 3);
            c_mn = val(a, i, 4);
        }
    }

    result = (-a_mn + b_mn*T + c_mn*pow(T,2))/T;
    return result;
}

/**
 * \f[
 * \Theta_m = \frac{Q_m X_m}{\sum_n Q_n X_n}
 * \f]
 * @param m Group ID
 * @param s Solution
 * @returns Theta_m
 */
double _Theta(int m, unifac_solution *s)
{
    int n;
    double sum = 0, result;
    for(n=0; n<s->dat->ngroups; n++)
        sum += s->dat->rows[n].Qk * _Xm(n, s);
    result = s->dat->rows[m].Qk * _Xm(m, s) / sum;
    return result;
}


/**
 * Activity of an isolated group of ID k in solution s.
 * \f[
 * \ln\Gamma_k = Q_k\left[1 - \ln\sum_m\Theta_m\Psi_{mk}
 *      - \sum_m\frac{\Theta_m\Psi_{km}}{\sum_n\Theta_n\Psi_{nm}}\right]
 * \f]
 * @param k Group ID
 * @param s Solution
 * @returns Gamma_k
 */
double _ln_Gamma(int k, unifac_solution *s, double T)
{
    int m, n;
    double Qk = s->m[k].dat->Qk,
           sum1 = 0, sum2 = 0, sum3 = 0,
           result;

    for(m=0; m<s->dat->ngroups; m++)
        sum1 += _Theta(m, s) * _Psi(m, k, s, T);

    for(m=0; m<s->dat->ngroups; m++) {
        for(n=0; n<s->dat->ngroups; n++)
            sum3 += _Theta(n, s) * _Psi(n, m, s, T);
        sum2 += _Theta(m, s)*_Psi(k, m, s, T)/sum3;
        sum3 = 0;
    }

    result = Qk * (1 - log(sum1) - sum2);
    return result;
}

double _ln_GammaElec(int k, unifac_solution *s, double T)
{
    int m, n;
    double Qk,
           sum1 = 0, sum2 = 0, sum3 = 0,
           result;

    for(m=0; m<s->dat->ngroups; m++) {
        sum1 += _Theta(m, s) * _PsiElec(m, k, s, T);
    }

    for(m=0; m<s->dat->ngroups; m++) {
        for(n=0; n<s->dat->ngroups; n++)
            sum3 += _Theta(n, s) * _PsiElec(n, m, s, T);

        if(sum3 != 0)
            sum2 += _Theta(m, s)*_PsiElec(k, m, s, T)/sum3;

        sum3 = 0;
    }
    if(solution_group_count(k, s))
        Qk = s->dat->rows[k].Qk;
    else
        Qk = 0;

    if(sum1 == 0)
        return 0;
    else
        return Qk * (1 - log(sum1) - sum2);

    //return result;
}

double _ln_gammar(int i, unifac_solution *s, double T)
{
    int k;
    double sum = 0,
           lnGammak, lnGammaki;
    unifac_solution *spure;
    for(k=0; k<s->dat->ngroups; k++) {
        lnGammak = _ln_Gamma(k, s, T);
        spure = UnifacPureSolution(i, s);
        lnGammaki = _ln_Gamma(k, s, T);
        UnifacDestroySolution(spure);
        sum += molec_group_count(k, s->m[i]) * (lnGammak - lnGammaki);
    }
    return sum;
}

double _ln_gammarElec(int i, unifac_solution *s, double T)
{
    int k;
    double sum = 0,
           lnGammak, lnGammaki;
    unifac_solution *spure;
    for(k=0; k<s->dat->ngroups; k++) {
        lnGammak = _ln_GammaElec(k, s, T);
        spure = UnifacPureSolution(i, s);
        lnGammaki = _ln_GammaElec(k, spure, T); 
        UnifacDestroySolution(spure);
        sum += molec_group_count(k, s->m[i]) * (lnGammak - lnGammaki);
    }
    return sum;
}

double ln_gamma(int i, unifac_solution *s, double T)
{
    if(nCols(s->dat->interactions) == 3)
        return _ln_gammac(i, s) + _ln_gammar(i, s, T);
    else
        return _ln_gammacElec(i, s) + _ln_gammarElec(i, s, T);
}

