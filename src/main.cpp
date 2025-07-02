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

    /***************************** Define PT Parameters *****************************/
    const auto vw = 0.7;
    const auto alpha = 0;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    const auto wN = PhaseTransition::dflt_PTParams::wN;
    const auto model = PhaseTransition::dflt_PTParams::model;
    // const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    PhaseTransition::Universe un;
    PhaseTransition::PTParams params1(vw, alpha, beta, dtau, wN, model, "exp", un);
    PhaseTransition::PTParams params2(vw, alpha, beta, dtau, wN, model, "sim", un);
    
    params1.print();
    /********************************************************************************/

    // Ex 1: Construct fluid profile
    Hydrodynamics::FluidProfile prof(params1);
    prof.write("test_prof.csv");
    prof.plot("test_prof.png");

    // Ex 2: Construct GW Spectrum
    const auto kRs_vec = logspace(1e-3, 1e+3, 5);
    const auto OmegaGW1 = Spectrum::GWSpec(kRs_vec, params1);
    // const auto OmegaGW2 = Spectrum::GWSpec(kRs_vec, params2);
    
    // add .plot() function for GWSpec!!
    plt::figure_size(800, 600);
    plt::loglog(OmegaGW1.kvec(), OmegaGW1.Pvec(), "k-");
    // plt::loglog(OmegaGW2.kvec(), OmegaGW2.Pvec(), "r-");
    plt::xlabel("K=kRs");
    plt::ylabel("Omega_GW(K)");
    plt::xlim(1e-3, 1e+3);
    // plt::ylim(1e-4, 1e+0);
    plt::grid(true);
    plt::save("GW_spectrum.png");

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    ProfilerStop();
    /******************************************************************/
    return 0;
}
