#ifndef __NORMALIZATION_H__
#define __NORMALIZATION_H__

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

int_t pow(int_t base, int_t exp);

class Normalization
{
public:
    typedef ring_int_real<int_t> RootTwoRing;
    typedef ring_int<int_t> OmegaRing;

    Normalization(const OmegaRing &z_poly);
    matrix2x2<int_t> get_result() const;
    void solve();

private:
    OmegaRing solve_y();
    RootTwoRing build_norm_equation();
    bool find_valid_r();
    real_t RootTwoToReal(const RootTwoRing &r);

    inline bool in_range(const RootTwoRing &r)
    {
        return RootTwoToReal(r) >= x_min && RootTwoToReal(r) <= x_max && abs(RootTwoToReal(r.g_conjugate())) <= r_dot_upper_bound;
    }

    inline bool simple_check_solvable(const RootTwoRing &rz)
    {
        return RootTwoToReal(rz) >= 0 && RootTwoToReal(rz.g_conjugate()) >= 0;
    }

    inline real_t max(const real_t &a, const real_t &b) { return a > b ? a : b; }

    OmegaRing z_poly;
    int_t L1;
    real_t ita;
    int_t R0;
    real_t r_dot_upper_bound;
    real_t x_min;
    real_t x_max;
    RootTwoRing r;
    int_t Lr;
    OmegaRing res_y;
    OmegaRing res_z;
};

#endif // __NORMALIZATION_H__