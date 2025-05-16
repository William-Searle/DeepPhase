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

    const PhaseTransition::PTParams params;
    // const auto kRs_vals = logspace(1e-3, 1e+3, 5);
    // const auto OmegaGW = Spectrum::GWSpec(kRs_vals, params);
    // const auto Ek = Spectrum::Ekin(kRs_vals, params);
    
    // const bool plot = false;
    // const bool plot = true;
    // if (plot) {
    //     plt::figure_size(800, 600);
    //     plt::loglog(OmegaGW.kvec(), OmegaGW.Pvec());
    //     plt::xlabel("K=kRs");
    //     plt::ylabel("Omega_GW(K)");
    //     // plt::xlim(1e-1, 1e+3);
    //     // plt::ylim(1e-4, 1e+0);
    //     plt::grid(true);
    //     plt::save("GW_spectrum.png");
    // }

    const auto k_vals = logspace(1e-3, 1e+3, 5);
    const auto p_vals = linspace(1e-2, 1e+3, 200);
    const auto z_vals = linspace(-1.0, 1.0, 200);
    const auto dlta = Spectrum::dlt2(k_vals, p_vals, z_vals, params);

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    ProfilerStop();
    /******************************************************************/
    return 0;
}
