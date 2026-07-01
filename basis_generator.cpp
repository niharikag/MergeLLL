/*
* Copyright (C) 2026 Niharika
 * Licensed under the GNU General Public License v3.0 or later.
 * https://www.gnu.org/licenses/gpl-3.0.html
 */

#include <iostream>
#include <vector>
#include <random>
#include "common.h"


/**
* Generates NTRU-style lattice
*
* [ qI   0 ]
* [ H    I ]
*
* H is circulant from a random polynomial h
 */

MatInt  generate_ntru_lattice(int N, long long q,
                    std::mt19937_64& rng) {

    int dim = 2 * N;

    MatInt Basis(dim, VecInt(dim, 0));

    for (int i = 0; i < N; i++) {
        Basis[i][i] = TO_MPZ(q);
       // mpz_set_si(Basis[i][i].get_data(), q);
    }
    std::uniform_int_distribution<long long> dist(0, q - 1);

    std::vector<long long> h(N);

    for (int i = 0; i < N; i++)
        h[i] = dist(rng);

    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            Basis[N + r][c] = TO_MPZ(
                h[(c - r + N) % N]);
        }
    }

    for (int i = 0; i < N; i++)
        Basis[N + i][N + i] = 1;

    return Basis;
}


MatInt generate_subsetsum_basis(int n, int bits, std::mt19937_64& rng) {
    VecInt a(n);
    VecInt x(n);

    // ---- instance generation ----
    std::uniform_int_distribution<long long> dist_a(0, (1LL << bits) - 1);

    for (int i = 0; i < n; i++) {
        a[i] = static_cast<long>(dist_a(rng) );
    }
    std::uniform_int_distribution<int> dist_x(0, 1);
    for (int i = 0; i < n; i++) {
        x[i] = dist_x(rng);
    }

    Real S = 0;
    for (int i = 0; i < n; i++) {
        S += a[i] * x[i];
    }

    // ---- basis construction ----
    MatInt B(n + 1, VecInt(n + 2, 0));

    for (int i = 0; i < n; i++) {
        B[i][i] = 2;
        B[i][n] = n * a[i];
        B[i][n + 1] = 0;
    }

    for (int i = 0; i < n; i++) {
        B[n][i] = 1;
    }

    B[n][n] = n * S;
    B[n][n + 1] = 1;

    return B;
}