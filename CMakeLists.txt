cmake_minimum_required(VERSION 4.1)
project(MergeLLL)

set(CMAKE_CXX_STANDARD 20)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FPLLL REQUIRED IMPORTED_TARGET fplll)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)
pkg_check_modules(MPFR REQUIRED IMPORTED_TARGET mpfr)
pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)
find_package(Threads REQUIRED)

add_executable(MergeLLL main.cpp
        utils.cpp
        lll.cpp
        potlll.cpp
        potlll_merge.cpp
        mergelll.cpp
        blockwise_kz.cpp
        basis_generator.cpp
        basis_quality_check.cpp
        run_experiment.cpp
)

target_link_libraries(MergeLLL PRIVATE PkgConfig::FPLLL PkgConfig::GMP PkgConfig::MPFR PkgConfig::GMPXX Threads::Threads)
