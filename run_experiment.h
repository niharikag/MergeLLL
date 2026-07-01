//
// Created by Niharika Gauraha on 2026-06-30.
//

#ifndef MERGELLL_RUN_EXPERIMENT_H
#define MERGELLL_RUN_EXPERIMENT_H

#include "basis_generator.h"
#include "basis_quality_check.h"
#include "lll.h"
#include "potlll_merge.h"

struct ExperimentInputParams {
    int n;
    int bits;
    int max_size;
    int block_size;
    bool run_all;
    double delta;
    std::string_view gen_name;
};

struct ExperimentResult {
    int n;

    int lll_swaps;
    double lll_rhf;
    Real lll_oe;
    Real lll_time;

    int merge_swaps;
    int merge_total_swaps;
    double merge_rhf;
    Real merge_oe;
    Real merge_time;

    double bkz_rhf;
    double bkz_oe;
    double bkz_time;
};

void run_experiment(const ExperimentInputParams &input_params);
ExperimentResult run_individual_experiment(const ExperimentInputParams &input_params);

#endif //MERGELLL_RUN_EXPERIMENT_H
