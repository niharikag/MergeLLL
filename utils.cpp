/*
* Copyright (C) 2026 Niharika
 * Licensed under the GNU General Public License v3.0 or later.
 * https://www.gnu.org/licenses/gpl-3.0.html
 */

# include "utils.h"
#include <gmpxx.h>

void print_matrix(const MatInt& A)
{
    for (const auto& row : A)
    {
        for (const auto& x : row)
        {
            std::cout << x << " ";
        }
        std::cout << "\n";
    }
}

// dot product
Real dot(const VecReal &a, const VecReal &b) {
    Real s = 0;
    for (int i = 0; i < a.size(); i++)
        s += a[i] * b[i];
    return s;
}

// the dot product is converted to Real to avoid overflow
Real dot(const VecInt &a, const VecInt &b) {
    Real s = 0;
    for (int i = 0; i < a.size(); i++) {
        Real ai = a[i];
        Real bi = b[i];
        s += ai * bi;
    }
    return s;
}

Real dot(const VecInt& a, const VecReal& b) {
    Real res = 0;
    for (int i = 0; i < a.size(); ++i) {
        res += static_cast<Real>(a[i]) * b[i];
    }
    return res;
}

// vector subtraction: a = a - factor * b
VecInt sub(VecInt a, const VecInt &b, int factor) {
    for (int i = 0; i < a.size(); i++)
        a[i] -= factor * b[i];
    return a;
}

VecReal sub(VecReal a, const VecReal &b, Real factor) {
    for (int i = 0; i < a.size(); i++)
        a[i] -= factor * b[i];
    return a;
}


void gram_schmidt(
    const MatInt& Basis,
    MatReal& Bstar,
    MatReal& mu,
    VecReal& B) {
    const int n = Basis.size();

    mu.assign(n, std::vector<Real>(n, 0));

    for (auto i = 0; i < n; i++) {
        Bstar[i] = VecReal(Basis[i].begin(), Basis[i].end());

        for (int j = 0; j < i; j++) {
             // if (B[j] < EPS) {
             //     B[j] = EPS;
             // }
            mu[i][j] = dot(Basis[i], Bstar[j]) / B[j];
            Bstar[i] = sub(Bstar[i], Bstar[j], mu[i][j]);
        }
        B[i] = dot(Bstar[i], Bstar[i]);
    }
    for (int i = 0; i < B.size(); ++i) {
        if (B[i] == 0)
            std::cerr << "B[" << i << "] is zero\n";
    }
}
mpz_class round_mpf(const mpf_class& x) ;

void size_reduce(MatInt &Basis, MatReal &mu, const int k, const int l) {
    //Int r = static_cast<int> (std::round(mu[k][l]));
    //int r = static_cast<int>(std::llround( mu[k][l].get_d()));
    //long long r = std::llround( mu[k][l].get_d());

    mpz_class r(round_mpf(mu[k][l]));

    if (r == 0)
        return;

    //Basis[k] -= r * Basis[l];
    auto &bk = Basis[k];
    const auto &bl = Basis[l];

    for (int j = 0; j < bk.size(); j++) {
        bk[j] -= __gmp_expr<mpz_t, mpz_t>(r) * bl[j];
    }

    mu[k][l] -= r;

    for (int j = 0; j <l; j++) {
        mu[k][j] -=  r * mu[l][j];

    }

}


bool lovasz(const MatReal& mu, const VecReal& B, const int k, const double delta) {
    return B[k] >= (delta - mu[k][k - 1] * mu[k][k - 1]) * B[k - 1];
}

#include <gmpxx.h>

mpz_class round_mpf(const mpf_class& x) {
    mpf_class fx = floor(x);
    mpf_class diff = x - fx;

    if (diff >= 0.5)
        return mpz_class(fx + 1);

    return mpz_class(fx);
}