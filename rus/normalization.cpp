#include <numeric>
#include <iostream>
#include <gmpxx.h>
#include <complex>
#include <random>
#include <cassert>
#include <map>
#include <mpfr.h>

#include "numeric_definition.h"
#include "rint.h"
#include "../appr/normsolver.h"
#include "../es/exactdecomposer.h"
#include "normalization.h"

namespace
{
real_t mpf_log2(const real_t &x)
{
    mpfr_t xx;
    mpfr_t yy;
    mpfr_init(xx);
    mpfr_init(yy);
    mpfr_set_f(xx, x.get_mpf_t(), MPFR_RNDN);
    mpfr_log2(yy, xx, MPFR_RNDN);
    real_t result;
    mpfr_get_f(result.get_mpf_t(), yy, MPFR_RNDN);
    mpfr_clear(xx);
    mpfr_clear(yy);
    return result;
}

int_t mpf_ceil_int(const real_t &x)
{
    real_t c;
    mpf_ceil(c.get_mpf_t(), x.get_mpf_t());
    return int_t(c);
}

real_t root_two_to_real(const Normalization::RootTwoRing &r)
{
    return r[0] + r[1] * SQRT2;
}

struct NormInitParams
{
    int_t L1;
    real_t ita;
    int_t R0;
    real_t r_dot_upper_bound;
};

NormInitParams compute_norm_init_params(const Normalization::OmegaRing &z_poly)
{
    const real_t z_abs2 = root_two_to_real(z_poly.abs2());
    const int_t L1 = mpf_ceil_int(mpf_log2(z_abs2));
    const real_t ita = real_t(L1) - mpf_log2(z_abs2);

    const real_t ratio = root_two_to_real(z_poly.g_conjugate().abs2()) / z_abs2;
    const real_t inner = ratio * real_t(L1) * real_t(L1) * NU * NU;
    const int_t R0 = mpf_ceil_int(mpf_log2(inner) / 2);

    const real_t exp_val = (real_t(R0) + ita) / 2;
    const int_t exp_int = int_t(exp_val);
    const real_t r_dot_upper_bound(real_t(pow(int_t(2), exp_int)));

    return {L1, ita, R0, r_dot_upper_bound};
}
} // namespace

int_t pow(int_t base, int_t exp)
{
    int_t result = 1;
    while (exp)
    {
        if (exp % 2 == int_t(1))
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

Normalization::Normalization(const OmegaRing &z_poly)
    : z_poly(z_poly),
      L1(0),
      ita(0),
      R0(0),
      r_dot_upper_bound(0),
      x_min(0),
      x_max(0),
      r(RootTwoRing(0, 0)),
      Lr(0)
{
    const NormInitParams p = compute_norm_init_params(z_poly);
    L1 = p.L1;
    ita = p.ita;
    R0 = p.R0;
    r_dot_upper_bound = p.r_dot_upper_bound;
    x_min = (1 - 1 / (2 * real_t(L1))) * r_dot_upper_bound;
    x_max = r_dot_upper_bound;
}

matrix2x2<int_t> Normalization::get_result() const
{
    assert(Lr >= 0);
    return matrix2x2<int_t>(
        res_z, res_y,
        -res_y.conjugate(), res_z.conjugate(),
        Lr.get_si());
}
void Normalization::solve()
{
    res_y = solve_y();
    res_z = z_poly * r;
}

Normalization::OmegaRing Normalization::solve_y()
{
    using std::abs;
    int_t max_try = 20;
    const normSolver &ns = normSolver::instance();
    OmegaRing y;
    while (--max_try >= 0)
    {
        RootTwoRing rz = build_norm_equation();
        // std::cout << "y^2 should be: " << rz << std::endl;
        // std::cout << "r is" << r << std::endl;

        if (!simple_check_solvable(rz))
            continue;

        if (ns.solve(rz, y))
        {
            return y;
        }
    }
    std::cerr << "Diverge" << std::endl;
    throw std::runtime_error("Diverge");
}
Normalization::RootTwoRing Normalization::build_norm_equation()
{
    int max_trial = 10;
    // r may not exist for Delta*delta = 2*r_dot_upper_bound * (x_max - x_min) < (1+\sqrt{2})^2
    while (!find_valid_r())
    {
        if (--max_trial < 0)
            throw std::runtime_error("No valid r found");
        x_min *= 2;
        x_max *= 2;
    }
    // std::cout << "r = " << r[0] << "+" << r[1] << "sqrt(2)" << std::endl;
    x_min *= 2;
    x_max *= 2;

    RootTwoRing rz = z_poly.abs2() * r.abs2();

    const real_t abs_sqr = RootTwoToReal(rz);
    Lr = mpf_ceil_int(mpf_log2(abs_sqr));
    rz = RootTwoRing(pow(2, Lr), 0) - rz;
    return rz;
}

// todo: r seems not so strict? Maybe we can use a better way to find r
bool Normalization::find_valid_r()
{
    // const real_t Delta = 2 * r_dot_upper_bound;
    real_t a_real(real_t(mpf_ceil_int((x_min + r_dot_upper_bound) / 2)));
    real_t b_real(real_t(mpf_ceil_int((x_min - r_dot_upper_bound) / (2 * SQRT2))));
    assert(a_real.fits_sint_p());
    assert(b_real.fits_sint_p());
    int_t a = a_real.get_si();
    int_t b = b_real.get_si();

    r = RootTwoRing(a, b);

    if (in_range(r))
        return true;

    ++b;
    r = RootTwoRing(a, b);
    if (in_range(r))
        return true;

    --b;
    --a;

    r = RootTwoRing(a, b);
    if (in_range(r))
        return true;

    return false;
}

real_t Normalization::RootTwoToReal(const RootTwoRing &r)
{
    return root_two_to_real(r);
}
