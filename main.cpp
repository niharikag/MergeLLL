#include <iostream>
#include <fplll.h>

#include "basis_generator.h"
#include "basis_quality_check.h"
#include "lll.h"
#include "potlll.h"
#include "potlll_merge.h"
#include "mergelll.h"
#include "run_experiment.h"
#include "blockwise_kz.h"
#include <fplll/nr/nr.h>

//using namespace fplll;


int run_individual_experiment(int n) ;

#include <fplll.h>

using namespace fplll;


int main() {
    FP_NR<mpfr_t>::set_prec(512);
    ExperimentInputParams input_params;
    input_params.gen_name="subsetsum";
    input_params.delta = DEFAULT_DELTA;
    input_params.n = 100;
    input_params.max_size=100;
    input_params.block_size=30;
    input_params.bits=20;
    //run_experiment(input_params);
    run_individual_experiment(input_params);
    return 0;
}

