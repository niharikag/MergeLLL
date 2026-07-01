/*
* Copyright (C) 2026 Niharika
 * Licensed under the GNU General Public License v3.0 or later.
 * https://www.gnu.org/licenses/gpl-3.0.html
 */

//#include <cmath>
#include <iostream>
# include "basis_quality_check.h"
#include "utils.h"

#include <gmpxx.h>
#include <mpfr.h>

void mpf_to_mpfr(const mpf_class& in, mpfr_t out, mpfr_rnd_t rnd = MPFR_RNDN) {
    // extract underlying GMP type
    mpf_srcptr tmp = in.get_mpf_t();

    // direct conversion (MPFR supports GMP float import)
    mpfr_set_f(out, tmp, rnd);
}


// check delta-LLL reduced
bool is_LLL_reduced(const MatInt& Basis, double delta) {
    auto n = Basis.size();
    auto m = Basis[0].size();

    MatReal Bstar(n, VecReal(m));
    MatReal mu(n, VecReal(n, 0.0));
    VecReal B(n, 0.0);

    gram_schmidt(Basis, Bstar, mu, B);


    bool size_ok = true;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            if (mu[i][j].get_d() > 0.5 + 1e-12) {
                std::cout << "Size FAIL mu[" << i << "," << j << "] = " << mu[i][j] << "\n";
                size_ok = false;
            }
        }
    }

    bool lovasz_ok = true;
    for (int k = 1; k < n; k++) {
        Real lhs = delta * B[k - 1];
        Real rhs = B[k] + mu[k][k - 1] * mu[k][k - 1] * B[k - 1];

        if (lhs > rhs) {
            std::cout << "Lovasz FAIL at k=" << k << "\n";
            lovasz_ok = false;
        }
    }

    std::cout << "Size reduced: " << std::boolalpha << size_ok << "\n";
    std::cout << "Lovasz condition: " <<std::boolalpha  << lovasz_ok << "\n";
    std::cout << "LLL reduced: " << std::boolalpha << (size_ok && lovasz_ok) << "\n";

    return size_ok && lovasz_ok;
}


// find shortest vector
VecInt shortest_vector_approx(const MatInt& Basis) {
    Real best_norm;
    //mpf_set_default_prec(512);
    bool has_best = false;
    VecInt best_vec;

    for (const auto& b : Basis) {
        Real norm_sq = dot(b, b);
        if (!has_best) {
            has_best = true;
            best_norm = norm_sq;
        }
        // std::cout << "norm_sq = " << norm_sq<< "\n";
        // std::cout << "best_norm = " << best_norm<< "\n";
        if (norm_sq <= best_norm) {
            best_norm = norm_sq;
            best_vec = b;
        }
    }

    return best_vec;
}


// compute log lattice volume
Real log_lattice_volume(const MatInt& Basis) {
    auto n = Basis.size();
    auto m = Basis[0].size();

    MatReal Bstar(n, VecReal(m));
    MatReal mu(n, VecReal(n, 0.0));
    VecReal B(n, 0.0);

    gram_schmidt(Basis, Bstar, mu, B);

    Real sum = 0.0;
    for (int i = 0; i < n; i++) {
        // if (B[i] > EPS)
        sum += std::log(std::sqrt(B[i].get_d()) );
    }
    return sum;
}


double root_hermite_factor(const MatInt& Basis) {
    const double n =Basis.size();

    Real log_vol = log_lattice_volume(Basis);
    //std::cout << "log_vol = " << log_vol<< "\n";
    VecInt best_vec = shortest_vector_approx(Basis);
    Real norm_sq = dot(best_vec, best_vec);
    //std::cout << "norm_sq = " << norm_sq<< "\n";
    double log_hf = 0.5L * std::log(norm_sq.get_d()) - log_vol.get_d() /n ;

    double rhf = std::exp(log_hf / n);

    //std::cout << "rhf = " << rhf<< "\n";
    return rhf;
}