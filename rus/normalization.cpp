#include <numeric>
#include <iostream>
#include <gmpxx.h>
#include <complex>
#include <random>
#include <cassert>
#include <map>

#include "numeric_definition.h"
#include "rint.h"
#include "../appr/normsolver.h"
#include "../es/exactdecomposer.h"
#include "normalization.h"

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

Normalization::Normalization(const OmegaRing &z_poly) : z_poly(z_poly),
                                                        L1(ceil(log2(RootTwoToReal(z_poly.abs2()).get_d()))),
                                                        ita(L1 - log2(RootTwoToReal(z_poly.abs2()).get_d())),
                                                        R0(ceil(
                                                            log2(
                                                                real_t(real_t(
                                                                           RootTwoToReal(z_poly.g_conjugate().abs2()) / RootTwoToReal(z_poly.abs2())) *
                                                                       L1 * L1 * NU * NU)
                                                                    .get_d()) /
                                                            2)),
                                                        r_dot_upper_bound(pow(2, real_t((R0 + ita) / 2).get_d())),
                                                        x_min((1 - 1 / (2 * real_t(L1))) * r_dot_upper_bound),
                                                        x_max(r_dot_upper_bound),
                                                        r(RootTwoRing(0, 0)),
                                                        Lr(0)
{
}
matrix2x2<int_t> Normalization::get_result() const
{
    assert(Lr >= 0);
    return matrix2x2<int_t>(
        res_z, res_y,
        -res_y.conjugate(), res_z.conjugate(),
        Lr.get_d());
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

    // assume it's small enough
    const double abs_sqr = RootTwoToReal(rz).get_d();
    Lr = ceil(log2(abs_sqr));
    rz = RootTwoRing(pow(2, Lr), 0) - rz;
    return rz;
}

// todo: r seems not so strict? Maybe we can use a better way to find r
bool Normalization::find_valid_r()
{
    // const real_t Delta = 2 * r_dot_upper_bound;
    real_t a_real = real_t(ceil((x_min + r_dot_upper_bound) / 2));
    real_t b_real = real_t(ceil((x_min - r_dot_upper_bound) / (2 * SQRT2)));
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
    return r[0] + r[1] * SQRT2;
}