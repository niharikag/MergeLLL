/*
* Copyright (C) 2026 Niharika
 * Licensed under the GNU General Public License v3.0 or later.
 * https://www.gnu.org/licenses/gpl-3.0.html
 */
/* MergeLLL: potLLL is used when merging */

#include "potlll_merge.h"
#include <future>


int merge_potlll_reduce(MatInt &Basis_orig,
    int left,
    const int mid,
    const int right,
    const double delta) {

    MatInt Basis(Basis_orig.begin() + left,
           Basis_orig.begin() + right);

    int k = mid - left +1;

    int swap_count = potlll_reduce(Basis, k, delta);

    // writ it back in the original matrix
    for (int i = 0; i < Basis.size(); i++) {
        Basis_orig[left + i] = Basis[i];
    }

    return swap_count;
}

int potlll_merge_reduce( MatInt& Basis, const double delta) {
    int total_swaps = 0;
    int n = Basis.size();

    int size = 1;
    int n_swaps = 0;

    while (size < n) {

        int left = 0;

        while (left < n - 1) {

            int mid = std::min(left + size - 1, n - 1);
            int right = std::min(left + 2 * size - 1, n - 1);

            if (mid < right) {
                n_swaps = merge_potlll_reduce(Basis, left, mid, right + 1, delta);
            }

            left += 2 * size;
        }

        size *= 2;
    }
    return n_swaps;
}


int parallel_potlll_merge(MatInt& Basis, int & total_swaps, double delta) {
    int n = Basis.size();
    total_swaps = 0;
    int level_swaps = 0;

    for (int size = 1; size < n; size *= 2) {

        std::vector<std::future<int>> tasks;

        for (int left = 0; left < n - 1; left += 2 * size) {

            int mid = std::min(left + size - 1, n - 1);
            int right = std::min(left + 2 * size - 1, n - 1);

            if (mid < right) {
                tasks.push_back(std::async(std::launch::async,
                    [&Basis, left, mid, right, delta]() -> int {
                        return merge_potlll_reduce(Basis, left, mid, right + 1, delta);
                    }
                ));
            }
        }

        //BARRIER: wait for all tasks in this level
        level_swaps = 0;
        for (auto& t : tasks) {
            level_swaps += t.get();
        }

        total_swaps += level_swaps;
    }

    return level_swaps;
}