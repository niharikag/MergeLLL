/*
* Copyright (C) 2026 Niharika
 * Licensed under the GNU General Public License v3.0 or later.
 * https://www.gnu.org/licenses/gpl-3.0.html
 */
/* MergeLLL: blockwiseKZ is used for reduction in each block then potLLL is used for merging */

#include "potlll_merge.h"
#include <future>
#include <fplll.h>
#include "blockwise_kz.h"

using namespace fplll;

void bkz_insert(MatInt &Basis, VecInt c);

int merge_using_potlll_reduce(MatInt &Basis_orig,
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


int mergelll_reduce( MatInt& Basis, int & total_swaps, const int block_size, const double delta) {
    total_swaps = 0;
    int n = Basis.size();

    int size = block_size;
    int n_swaps = 0;

    for (int b = 0; b < n; b += block_size) {
        int end = std::min(b + block_size, n);
        MatInt subBasis(Basis.begin() + b,
           Basis.begin() + end);

        int status = bkz_wrapper(subBasis, subBasis.size(), delta);
        if (status==0) {
            // writ it back in the original matrix
            for (int i = 0; i < subBasis.size(); i++) {
                Basis[b + i] = subBasis[i];
            }
        }
        // std::cout << "print subbasis = " << "\n";
        // print_matrix(subBasis);
    }

    while (size < n) {

        int left = 0;

        while (left < n - 1) {

            int mid = std::min(left + size - 1, n - 1);
            int right = std::min(left + 2 * size - 1, n - 1);

            if (mid < right) {
                n_swaps = merge_using_potlll_reduce(Basis, left, mid, right + 1, delta);
                total_swaps += n_swaps;
            }

            left += 2 * size;
        }

        size *= 2;
    }
    return n_swaps;
}

int parallel_mergelll_reduce(MatInt& Basis, int & total_swaps, const int block_size, double delta) {
    int n = Basis.size();
    total_swaps = 0;
    int level_swaps = 0;
    std::vector<std::future<void>> block_task;

    for (int b = 0; b < n; b += block_size) {
        block_task.push_back(std::async(std::launch::async,
        [&, b]() {
            int end = std::min(b + block_size, n);
            MatInt subBasis(Basis.begin() + b,
               Basis.begin() + end);

            int status = bkz_wrapper(subBasis, subBasis.size(), delta);

            if (status==0) {
                for (int i = 0; i < subBasis.size(); i++) {
                    Basis[b + i] = subBasis[i];
                }
            }
        }));
    }

    // Wait for all the blocks to finish.
    for (auto &f : block_task) {
        f.get();
    }

    // Now combine the blocks using potLLL
    for (int size = block_size; size < n; size *= 2) {

        std::vector<std::future<int>> tasks;

        for (int left = 0; left < n - 1; left += 2 * size) {

            int mid = std::min(left + size - 1, n - 1);
            int right = std::min(left + 2 * size - 1, n - 1);

            if (mid < right) {
                tasks.push_back(std::async(std::launch::async,
                    [&Basis, left, mid, right, delta]() -> int {
                        return merge_using_potlll_reduce(Basis, left, mid, right + 1, delta);
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
