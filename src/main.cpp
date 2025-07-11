// main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <gperftools/profiler.h>
#include <omp.h>

// modify include list when testing of program finished - currently includes everything
#include "hydrodynamics.hpp"
#include "PhaseTransition.hpp"
#include "spectrum.hpp"
#include "profile.hpp"
#include "physics.hpp"
#include "maths_ops.hpp"
#include "constants.hpp"
#include "tests.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/*
TO DO:
- use libraries?
- Input file where you can specify universe and pt constants
- Write clock
- test spectrum functions (Ekin, variants of delta, GW spectrum) behaviour in limits they discuss in paper
*/

int main() {
    /************************ CLOCK / PROFILER *************************/
    ProfilerStart("profile.out");
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    /*** targets for k, p, z, tau: ***/
    // const auto kRs_vals = logspace(1e-3, 1e+3, 100);
    // const auto pRs_vals = logspace(1e-2, 1e+3, 1000);
    // const auto Ttilde_vals = logspace(1e-2, 20, 1000);
    // const auto z_vals = linspace(-1.0, 1.0, 1000);
    // takes 2-4mins to run
    /****************************/
    
    // example_Kin_Spec("Ekin_spectrum");
    // example_GW_Spec("../GW_spectrum");

    // test_SiCi_spline();

    const auto vw = 0.5;
    const auto alN = 0.1;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    const auto wN = PhaseTransition::dflt_PTParams::wN;
    const auto model = PhaseTransition::dflt_PTParams::model;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const PhaseTransition::Universe un;
    const PhaseTransition::PTParams params(vw, alN, beta, dtau, wN, model, nuc_type, un);

    un.print();
    params.print();

    int n = 2;
    const auto k_vals = linspace(1e-3, 1e+3, n);
    const auto p_vals = linspace(1e-2, 1e+3, n);
    const auto z_vals = linspace(-1.0, 1.0, n);

    const auto dlt = Spectrum::dlt_SSM(k_vals, p_vals, z_vals, params);

    for (int kk = 0; kk < k_vals.size(); kk++) {
        const auto k = k_vals[kk];
        for (int pp = 0; pp < p_vals.size(); pp++) {
            const auto p = p_vals[pp];
            for (int zz = 0; zz < z_vals.size(); zz++) {
                const auto z = z_vals[zz];
                std::cout << "k,p,z=" << k_vals[kk] << "," << p_vals[pp] << "," << z_vals[zz] << "\n"
                          << "dlt=" << dlt[kk][pp][zz] << "\n\n";

            }
        }
    }
    

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    ProfilerStop();
    /******************************************************************/
    return 0;
}
