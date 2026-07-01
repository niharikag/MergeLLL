/*
* Copyright (C) 2026 Niharika
 * Licensed under the GNU General Public License v3.0 or later.
 * https://www.gnu.org/licenses/gpl-3.0.html
 */

#include "potlll.h"
#include <vector>
#include <cmath>
#include <limits>


static void deep_insert(MatInt& Basis, int j, int k) {
    VecInt tmp = std::move(Basis[k]);

    for (int i = k; i > j; --i)
        Basis[i] = std::move(Basis[i - 1]);

    Basis[j] = std::move(tmp);
}

static Real insertion_ratio(const MatReal& mu, const VecReal& B, const int j, const int k) {
    Real ratio = 1.0;
    auto proj = B[k];
    for (int i = 0; i <= k; ++i) {
        if (B[i] == 0)
            std::cerr << "In Insetion ration, B[" << i << "] is zero\n";
    }
    for (int l = k - 1; l >= j; --l) {
        proj += mu[k][l] * mu[k][l] * B[l];
        ratio *= proj / B[l];
    }
    if (ratio <=0)
        std::cerr << "ratio" << "\n";
    return ratio;
}
// compute Gram Schmidt for kth vector
static void compute_gs_row(MatInt& Basis, MatReal& Bstar, MatReal& mu, VecReal& B, const int k) {
    Bstar[k] = VecReal(Basis[k].begin(), Basis[k].end());
    //Bstar[k] = Basis[k];

    for (int j = 0; j < k; j++) {
         // if (B[j] < EPS) {
         //     B[j] = EPS;
         // }

        if (B[j] == 0) {
            std::cerr << "B[" << j << "] is zero\n";
            for (int s = 0; s < Basis[j].size(); s++) {
                std::cerr << Basis[j][s] << " "<< "\n";
            }
        }
        else {
            mu[k][j] = dot(Basis[k], Bstar[j]) / B[j];
        }
        Bstar[k] = sub(Bstar[k], Bstar[j], mu[k][j]);
    }

    B[k] = dot(Bstar[k], Bstar[k]);

}

//compute GS for first k vectors
void gram_schmidt_upto_k(
    const MatInt& Basis,
    MatReal& Bstar,
    MatReal& mu,
    VecReal& B,
    const int k) {

    for (auto i = 0; i <= k; i++) {
        Bstar[i] = VecReal(Basis[i].begin(), Basis[i].end());

        for (int j = 0; j < i; j++) {
            if (B[j] == 0)
                std::cerr << "B[" << j << "] is zero\n";
            mu[i][j] = dot(Basis[i], Bstar[j]) / B[j];
            Bstar[i] = sub(Bstar[i], Bstar[j], mu[i][j]);
        }
        B[i] = dot(Bstar[i], Bstar[i]);
        if (B[i] == 0) {
            for (int s = 0; s < Bstar[i].size(); s++) {
                std::cerr << Bstar[i][s] << " "<< "\n";
            }
            std::cerr << "B[" << i << "] is zero\n";
        }
    }
}

int potlll_reduce(MatInt &Basis, const int start_pos, const double delta) {

    int n = Basis.size();
    int m = Basis[0].size();

    int swap_count = 0;

    MatReal Bstar(n, VecReal(m));
    MatReal mu(n, VecReal(n, 0.0));
    VecReal B(n, 0.0);

    //gram_schmidt(Basis, Bstar, mu, B);

    int k = start_pos;
    gram_schmidt_upto_k(Basis, Bstar, mu, B, k);

    while (k < n) {

        compute_gs_row(Basis, Bstar, mu, B, k);

        for (int j = k - 1; j >= 0; --j) {
            size_reduce(Basis, mu, k, j);
        }

        //compute_gs_row(Basis, Bstar, mu, B, k);

        int best_j = k;
        Real best_ratio = 1.0;

        for (int j = 0; j < k; j++) {
            Real r = insertion_ratio(mu, B, j, k);

            if (r> 0 && r < best_ratio) {
                best_ratio = r;
                best_j = j;
            }
        }
        // std::cout << delta << " ";
        // std::cout << best_ratio << " " << std::endl;

        if (delta > best_ratio) {

            deep_insert(Basis, best_j, k);
            compute_gs_row(Basis, Bstar, mu, B, best_j);
            k = std::max(best_j, 1);

            swap_count++;

        } else {
            k++;
        }
    }

    return swap_count;
}
