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

/*
TO DO:
- use libraries?
- Input file where you can specify universe and pt constants
- Write clock
- test spectrum functions (Ekin, variants of delta, GW spectrum) behaviour in limits they discuss in paper
*/

int main() {
    ProfilerStart("profile.out");
    const auto ti = std::chrono::high_resolution_clock::now();

    // const PhaseTransition::PTParams params;

    /*** targets for k, p, z, tau: ***/
    // const auto k_vals = logspace(1e-3, 1e+3, 100);
    // const auto p_vals = logspace(1e-2, 1e+3, 1000);
    // const auto Ttilde_vals = logspace(1e-2, 20, 1000);
    // const auto z_vals = linspace(-1.0, 1.0, 4000);
    // takes 2-4mins to run
    /****************************/

    // const auto dlta = Spectrum::dlt(k_vals, p_vals, z_vals, params);

    test_GWSpec(false);

    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    ProfilerStop();
    return 0;
}
