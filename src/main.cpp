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

    // const PhaseTransition::Universe un;
    // un.print();

    // const PhaseTransition::PTParams params;
    // params.print();

    test_GWSpec();
    // test_Ekin();

    // const auto k_vals = logspace(1e-3, 1e+3, 5);
    // const auto p_vals = linspace(1e-2, 1e+3, 200);
    // const auto z_vals = linspace(-1.0, 1.0, 200);

    // const auto nk = k_vals.size();
    // const auto np = p_vals.size();
    // const auto nz = z_vals.size();

    // const auto Ek = Spectrum::Ekin(kRs_vals, params);
    // const auto OmegaGW = Spectrum::GWSpec(kRs_vals, params);
    // const auto dlta1 = Spectrum::dlt(50, k_vals, p_vals, z_vals, params);
    // const auto dlta2 = Spectrum::dlt(150, k_vals, p_vals, z_vals, params);

    // int count = 0;
    // for (int kk = 0; kk < nk; kk++) {
    //     const auto k = k_vals[kk];
    //     for (int pp = 0; pp < np; pp++) {
    //         const auto p = p_vals[pp];
    //         for (int zz = 0; zz < nz; zz++) {
    //             const auto z = z_vals[zz];
    //             const auto diff = (dlta1[kk][pp][zz] - dlta2[kk][pp][zz]) / dlta1[kk][pp][zz];
    //             std::cout << "(k,p,z)=(" << k << "," << p << "," << z << "), diff=" << diff << "\n";
    //             if (diff > 0.1)
    //                 count++;
    //         }
    //     }
    // }
    // std::cout << "count=" << count << "\n";

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    ProfilerStop();
    /******************************************************************/
    return 0;
}
