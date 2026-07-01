//
// Created by Niharika Gauraha on 2026-06-29.
//

#ifndef MERGELLL_POTLLL_MERGE_H
#define MERGELLL_POTLLL_MERGE_H

#include "utils.h"
#include "potlll.h"

int potlll_merge_reduce( MatInt& Basis, double delta);
int parallel_potlll_merge_reduce( MatInt& Basis, int & total_swaps, double delta);

#endif //MERGELLL_POTLLL_MERGE_H
