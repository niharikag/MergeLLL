//
// Created by Niharika Gauraha on 2026-06-30.
//

#ifndef MERGELLL_BASIS_GENERATOR_H
#define MERGELLL_BASIS_GENERATOR_H
#include "common.h"

MatInt generate_ntru_lattice(int N, long long q=2147483647, std::mt19937_64& rng=g_rng);
MatInt generate_subsetsum_basis(int n, int bits, std::mt19937_64& rng);
#endif //MERGELLL_BASIS_GENERATOR_H
