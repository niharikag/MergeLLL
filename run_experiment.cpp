/*
* Copyright (C) 2026 Niharika
 * Licensed under the GNU General Public License v3.0 or later.
 * https://www.gnu.org/licenses/gpl-3.0.html
 */

#include "run_experiment.h"
#include <iostream>
#include <chrono>

#include "blockwise_kz.h"
#include "mergelll.h"
#include "potlll.h"
#include "lll.h"

// run experiment for each individual dimension
ExperimentResult run_individual_experiment(const ExperimentInputParams& input_params) {
    MatInt Basis_orig;

    if (input_params.gen_name == "subsetsum") {
        Basis_orig = generate_subsetsum_basis(input_params.n, input_params.bits, g_rng);
    } else {
        Basis_orig = generate_ntru_lattice(input_params.n/2, 2147483647, g_rng);
    }

    ExperimentResult res;
    res.n = Basis_orig.size();
    lll_reduce(Basis_orig, 1,input_params.delta);

    int total_swaps = 0;
    int n_swaps;

    {
        MatInt Basis_lll = Basis_orig;
        auto start = std::chrono::high_resolution_clock::now();
        n_swaps = potlll_reduce(Basis_lll, 1, input_params.delta);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Time potLLL: " << elapsed.count() << " seconds\n"<<std::flush;

        res.lll_swaps = n_swaps;
        //is_LLL_reduced(Basis_lll, input_params.delta);
        res.lll_rhf = root_hermite_factor(Basis_lll);
        //res.lll_oe = compute_orthogonality_error(Basis_lll);
        res.lll_time = elapsed.count();
        std::cout << "potLLL output: " << n_swaps << " "<< res.lll_rhf << "\n"<<std::flush;
    }
    // MergeLLL
    {
        MatInt Basis_lll_merge = Basis_orig;
        auto start = std::chrono::high_resolution_clock::now();
        //print_matrix(Basis_lll_merge);
        n_swaps = parallel_mergelll_reduce(Basis_lll_merge, total_swaps, input_params.block_size, input_params.delta);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Time LLL Merge: " << elapsed.count() << " seconds\n"<<std::flush;

        //is_LLL_reduced(Basis_lll_merge, input_params.delta);

        res.merge_swaps = n_swaps;
        res.merge_total_swaps = total_swaps;
        res.merge_rhf = root_hermite_factor(Basis_lll_merge);
        //res.merge_oe = compute_orthogonality_error(Basis_lll_merge);
        res.merge_time = elapsed.count();
        std::cout << "LLL Merge output: " << n_swaps <<" " << total_swaps<< " "<< res.merge_rhf << "\n"<<std::flush;

        //print_matrix(Basis_lll_merge);
    }
    // BKZ
    {
        MatInt Basis_bkz = Basis_orig;
        auto start = std::chrono::high_resolution_clock::now();
        bkz_wrapper(Basis_bkz, input_params.block_size, input_params.delta);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Time LLL Merge: " << elapsed.count() << " seconds\n"<<std::flush;

        res.bkz_rhf = root_hermite_factor(Basis_bkz);
        //res.bkz_oe = compute_orthogonality_error(Basis_lll_merge);
        res.bkz_time = elapsed.count();
        std::cout << "BKZ output: "<< res.bkz_rhf << "\n"<<std::flush;
    }

    return res;
}

// run experiment for all sizes from 20 to max_size (and for a given bit weight for subsetsum)
void run_experiment(const ExperimentInputParams &input_params) {
    std::cout << input_params.gen_name << "\n" << std::flush;

    std::vector<ExperimentResult> results;


    for (int n = 20; n <= input_params.max_size; n += 20) {
        std::cout << n << "\n" << std::flush;
        ExperimentInputParams new_input_params = input_params;
        new_input_params.n = n;
        results.push_back(
            run_individual_experiment(new_input_params)
        );
    }

    // // print consolidated output
    // std::cout << "\nDim | LLL_swaps | LLL_HF | LLL_OrthErr | "
    //              "merge_swaps | merge_total | merge_HF | merge_OrthErr"
    //              "BKZ_hnf \n";
    //
    // for (const auto& r : results) {
    //     std::cout << r.n << " | "
    //               << r.lll_swaps << " | " << r.lll_rhf << " | " << r.lll_oe << " | "
    //               << r.merge_swaps << " | " << r.merge_total_swaps << " | " << r.merge_rhf << " | " << r.merge_oe
    //                 << " | " << r.bkz_rhf << " | " << r.bkz_oe << " | "
    //               << "\n";
    // }
    // print RHF separately
    std::cout << "\nDim | LLL_swaps | "
                 "MergeBKZ_swaps | MergeBKZ_total_swaps \n";

    for (const auto& r : results) {
        std::cout << r.n << " | "
                  << r.lll_swaps << " | "
                  << r.merge_swaps << " | "
                    << r.merge_total_swaps
                  << "\n";
    }

    // print RHF separately
    std::cout << "\nDim | LLL_RHF | "
                 "MergeBKZ_RHF |"
                "BKZ_RHF | \n";;

    for (const auto& r : results) {
        std::cout << r.n << " | "
                  << r.lll_rhf << " | "
                  << r.merge_rhf << " | "
                    << r.bkz_rhf
                  << "\n";
    }

    // print CPU time separately
    std::cout << "\nDim | LLL_CPU_Time | "
                "MergeBKZ_CPU_Time |"
                 "BKZ_CPU_Time | \n";

    for (const auto& r : results) {
        std::cout << r.n << " | "
                  << r.lll_time << " | "
                  << r.merge_time << " | "
                    << r.bkz_time
                  << "\n";
    }

}

