/*
* Copyright (C) 2026 Niharika
 * Licensed under the GNU General Public License v3.0 or later.
 * https://www.gnu.org/licenses/gpl-3.0.html
 */

#include "lll.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <ostream>


void static swap_step(MatInt& Basis, MatReal& mu, VecReal& B, int k) {
    Real old_mukkm1 = mu[k][k - 1];

    Real old_Bkm1 = B[k - 1];
    Real old_Bk = B[k];

    Real B_temp = old_Bk + old_mukkm1 * old_mukkm1 * old_Bkm1;

    std::swap(Basis[k - 1], Basis[k]);

    B[k - 1] = B_temp;
    B[k] = old_Bkm1 * old_Bk / B_temp;

    mu[k][k - 1] = old_mukkm1 * old_Bkm1 / B_temp;

    for (int j = 0; j < k - 1; ++j) {
        Real mukj = mu[k][j];
        Real mukm1j = mu[k - 1][j];

        mu[k - 1][j] = mukj;
        mu[k][j] = mukm1j;
    }

    for (int i = k + 1; i < B.size(); ++i) {
        Real mik = mu[i][k];
        Real mikm1 = mu[i][k - 1];

        mu[i][k - 1] = mu[k][k - 1] * mikm1 + mik - mik * old_mukkm1 * mu[k][k - 1];
        mu[i][k] = mikm1 - old_mukkm1 * mik;
    }
}


int lll_reduce(MatInt &Basis, int start_pos, double delta) {
    int n = Basis.size();
    int m = Basis[0].size();

    MatReal Bstar(n, VecReal(m));
    MatReal mu(n, VecReal(n, 0.0));
    VecReal B(n, 0.0);

    int n_swaps = 0;

    gram_schmidt(Basis, Bstar, mu, B);

    int k = start_pos;

    while (k < n) {
        size_reduce(Basis, mu, k, k-1);
        //std::cout << "k+1 = " << k+1 << std::endl;
        if (lovasz(mu, B, k, delta)) {
            for (int j = k - 2; j >= 0; --j) {
                size_reduce(Basis, mu, k, j);
            }
            k++;
        } else {
            swap_step(Basis, mu, B, k);
            k = std::max(k - 1, 1);
            n_swaps++;
        }
    }

    return n_swaps;
}


int test_lll() {
    MatInt Basis = {
        {19, 2, 32, 46, 3, 33},
        {15, 42, 11, 0, 3, 24},
        {43, 15, 0, 24, 4, 16},
        {20, 44, 44, 0, 18, 15},
        {0, 48, 35, 16, 31, 31},
        {48, 33, 32, 9, 1, 29}
    };

    lll_reduce(Basis);
    print_matrix(Basis);
    // mpz_class x;
    // std::cout << x;
    return 0;
}