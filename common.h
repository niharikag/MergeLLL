//
// Created by Niharika Gauraha on 2026-06-29.
//

#ifndef MERGELLL_COMMON_H
#define MERGELLL_COMMON_H
#include <random>
#include <gmpxx.h>
#include <iostream>
#include <fplll/nr/matrix.h>

using Int = mpz_class;
using VecInt = std::vector<Int>;
using MatInt = std::vector<VecInt>;
using Real = mpf_class;
using VecReal = std::vector<Real>;
using MatReal = std::vector<VecReal>;
// using MatInt = fplll::ZZ_mat<mpz_class>;
// using MatReal = fplll::ZZ_mat<Real>;

using Float = long double;



// global constants
inline std::mt19937_64 g_rng(21);
//inline const Real EPS = 1e-12;
inline constexpr double DEFAULT_DELTA = 0.99;
inline constexpr int DEFAULT_START_POS = 1;

#define TO_MPZ(q) mpz_class(static_cast<long>(q))

#endif //MERGELLL_COMMON_H
