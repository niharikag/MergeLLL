//
// Created by Niharika Gauraha on 2026-06-30.
//
#include <vector>
#include <gmpxx.h>
#include "blockwise_kz.h"
#include <fplll.h>


using namespace fplll;

void matint_to_zz_mat(const MatInt &in, ZZ_mat<mpz_t> &out)
{
    const int rows = static_cast<int>(in.size());
    const int cols = rows ? static_cast<int>(in[0].size()) : 0;

    out.resize(rows, cols);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            mpz_set(out[i][j].get_data(), in[i][j].get_mpz_t());
        }
    }
}

void zz_mat_to_matint(const ZZ_mat<mpz_t> &in, MatInt& out)
{
    const int rows = static_cast<int>(in.get_rows());
    const int cols = rows ? static_cast<int>(in.get_cols()) : 0;

    //out.resize(rows);
    for (int i = 0; i < rows; i++)
    {
        //out[i].resize(cols);
        for (int j = 0; j < cols; j++)
        {
            mpz_set(out[i][j].get_mpz_t(), in[i][j].get_data());
        }
    }
}


int bkz_wrapper(MatInt &Basis, const int block_size, const double delta)
{
    ZZ_mat<mpz_t> Basis_zz; // convert MatInt to ZZ Mat
    matint_to_zz_mat(Basis, Basis_zz);

    ZZ_mat<mpz_t> U;

    // Load strategies
    std::vector<Strategy> strategies =
        load_strategies_json(strategy_full_path(default_strategy()));

    // BKZ parameters
    BKZParam param(block_size, strategies, delta);

    // Run BKZ
    bkz_reduction(&Basis_zz, &U, param, fplll::FT_DEFAULT);

    zz_mat_to_matint(Basis_zz,Basis );

    return 0;
}