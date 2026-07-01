//
// Created by Niharika Gauraha on 2026-06-29.
//

#ifndef MERGELLL_MERGE_LLL_H
#define MERGELLL_MERGE_LLL_H

#include "utils.h"

int mergelll_reduce(MatInt &Basis,  int& total_swaps, int block_size = 5, double delta = DEFAULT_DELTA);
int parallel_mergelll_reduce(MatInt& Basis, int& total_swaps, int block_size = 5, double delta = DEFAULT_DELTA);

#endif //MERGELLL_MERGE_LLL_H
