/**
 * @file unifac.c
 * Implementation of the UNIFAC method for predicting VLE.
 */

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
double _ri(unifac_molec *molec)
{
    int k, sum = 0;
    for(k=0; k<molec.ngroups; k++)
        sum += molec->ngroups * molec->ngroups->dat->Rk;
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
double _qi(unifac_molec *molec)
{
    int k, sum = 0;
    for(k=0; k<molec.ngroups; k++)
        sum += molec->ngroups * molec->ngroups->dat->Qk;
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
double _Li(unifac_molec *molec)
{
    double ri, qi;
    ri = _ri(molec);
    qi = _qi(molec);
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
double _phii(int i, unifac_solution *s)
{
    int j;
    double numer, denom = 0;
    numer = s->xi[i] * _ri(s->m[i]);
    for(j=0; j<s->nsolutes; j++)
        denom += s->xi[j] * _ri(s->m[j]);
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
double _thetai(int i, unifac_solution *s)
{
    int j;
    double numer, denom = 0;
    numer = s->xi[i] * _qi(s->m[i]);
    for(j=0; j<s->nsolutes; j++)
        denom += s->xi[j] * _qi(s->m[j]);
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
double _ln_gammaic(int i, unifac_solution *s)
{
    int j;
    double result, sum = 0;

    for(j=0; j<m->nsolutes; j++)
        sum += s->xi[i] * _Li(s->m[i]);

    result = log(_phii(i, s)/s->xi[i])
            + Z/2*_qi(s->m[i])*log(_thetai(i, s)/_phii(i, s))
            + _Li(s->m[i])
            - _phii(i, s)/s->xi[i] * sum;

    return result;
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
    double numer = 0, denom = 0;
    for(j=0; j<s->nsolutes; j++)
        numer += molec_group_count(m, s->m[j]) * s->xi[j];
    for(j=0; j<s->nsolutes; j++) {
        for(n=0; n<s->dat->ngroups; n++)
            denom += molec_group_count(n, s->m[j]) * s->xi[j];
    }

    return numer/denom;
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
double _Psimn(int m, int n, unifac_solution *s, double T)
{
    int id1 = s->m[m]->dat->id,
        id2 = s->m[n]->dat->id;
    double a = val(s->dat->interactions, id1, id2);
    return -a/T;
}

/**
 * \f[
 * \Theta_m = \frac{Q_m X_m}{\sum_n Q_n X_n}
 * \f]
 * @param m Group ID
 * @param s Solution
 * @returns Theta_m
 */
double _Thetam(int m, unifac_solution *s)
{
    int n;
    double sum = 0;
    for(n=0; n<s->dat->ngroups; n++)
        sum += s->dat->rows[n]->Qk * _Xm(n);
    return s->dat->rows[n]->Qk * _Xm(n) / sum;
}

/**
 * Take a count of the number of groups with a particular ID in a specific
 * molecule. Return zero if the group isn't in there.
 * @param id Group ID
 * @param m Molecule
 * @returns Number of groups
 */
int molec_group_count(int id, unifac_molec *m)
{
    int i = 0;
    while(i<m->ngroups) {
        if(id == m->ids[i])
            return m->count[i];
        i++;
    }
    return 0;
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
double _ln_Gammak(int k, unifac_solution *s)
{
    int m, n;
    double Qk = s->m[k]->dat->Qk,
           sum1 = 0, sum2 = 0, sum3 = 0,
           result;

    for(m=0; m<s->dat->ngroups; m++)
        sum1 += _Thetam(m, s) * _Psimn(m, k, s);

    for(m=0; m<s->dat->ngroups; m++) {
        for(n=0; n<s->ngroups; n++)
            sum3 += _Thetam(n, s) * _Psimn(n, m, s);
        sum2 += _Thetam(m, s)*_Psimn(k, m, s)/sum3;
        sum3 = 0;
    }

    result = Qk * (1 - log(sum1) - sum2);
    return result;
}

double _ln_gammair(int i, unifac_solution *s)
{
    int k;
    double sum = 0,
           lnGammak, lnGammaki;
    for(k=0; k<s->dat->ngroups; k++) {
        lnGammak = _lnGammak(k, s);
        lnGammaki = 0; /* Make a new solution consisting only of molecule i */
        sum += molec_group_count(k, s->m[i]) * (lnGammak - lnGammaki);
    }
    return sum;
}

