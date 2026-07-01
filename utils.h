//
// Created by Niharika Gauraha on 2026-06-29.
//

#ifndef MERGELLL_UTILS_H
#define MERGELLL_UTILS_H


# include "common.h"

void print_matrix(const MatInt &B);

void gram_schmidt(
    const MatInt& Basis,
    MatReal& Bstar,
    MatReal& mu,
    VecReal& B
);

void size_reduce(MatInt &Basis, MatReal &mu, int k, int l);

bool lovasz(const MatReal& mu, const VecReal& B, int k, double delta);


VecReal sub(VecReal a, const VecReal& b, Real factor);
VecReal sub(VecInt a, const VecInt& b, Real factor);

Real dot(const VecInt& a, const VecInt& b); // long double to increase precision
Real dot(const VecReal& a, const VecReal& b);
Real dot(const VecInt& a, const VecReal& b);




#endif //MERGELLL_UTILS_H
