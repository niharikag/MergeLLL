//
// Created by Niharika Gauraha on 2026-06-30.
//

#ifndef MERGELLL_BASIS_QUALITY_CHECK_H
#define MERGELLL_BASIS_QUALITY_CHECK_H

#include "common.h"

bool is_LLL_reduced(const MatInt& Basis, double delta=DEFAULT_DELTA);

double root_hermite_factor(const MatInt& Basis);

#endif //MERGELLL_BASIS_QUALITY_CHECK_H
