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
    const auto alN = 0.1;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    const auto wN = PhaseTransition::dflt_PTParams::wN;
    const auto model = PhaseTransition::dflt_PTParams::model;
    // const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    PhaseTransition::Universe un;
    PhaseTransition::PTParams params1(vw, alN, beta, dtau, wN, model, "exp", un);
    // PhaseTransition::PTParams params2(vw, alN, beta, dtau, wN, model, "sim", un);
    
    un.print();
    params1.print();
    /********************************************************************************/

    // const auto k_vals = logspace(1e-3, 1e+3, 5);
    // const auto p_vals = logspace(1e-2, 1e+3, 200);
    // const auto z_vals = linspace(-1.0, 1.0, 200);

    // const auto dlt2 = Spectrum::dlt_SSM(k_vals, p_vals, z_vals, params1);

    // for (int kk = 0; kk < k_vals.size(); kk++) {
    //     std::cout << "dlt(" << k_vals[kk] << ")=" << dlt2[kk][0][0] << "\n";
    // }
    
    // Ex 1: Construct fluid profile
    Hydrodynamics::FluidProfile prof1(params1);
    // Hydrodynamics::FluidProfile prof2(params2);
    // prof.write("fluid_prof.csv");
    prof1.plot("fluid_prof.png");

    // Ex 2: Construct kinetic power spectrum
    // const auto kRs_vals = logspace(1e-1, 1e+3, 500);

    // const auto Ek1 = Spectrum::Ekin(kRs_vals, prof1);
    // const auto Eks1 = Spectrum::zetaKin(Ek1);

    // const auto Ek2 = Spectrum::Ekin(kRs_vals, prof2);
    // const auto Eks2 = Spectrum::zetaKin(Ek2);
    
    // plt::figure_size(800, 600);
    // plt::loglog(Eks1.k(), Eks1.P(), "k-");
    // plt::loglog(Eks2.k(), Eks2.P(), "r-");
    // plt::xlabel("K=kRs");
    // plt::ylabel("Ekin(K)");
    // plt::xlim(1e-1, 1e+3);
    // plt::ylim(1e-5, 1e+0);
    // plt::grid(true);
    // plt::save("../Ekin_spectrum.png");

    // Ex 3: Construct GW Spectrum
    // const auto kRs_vals = logspace(1e-3, 1e+3, 100);
    // const auto OmegaGW = Spectrum::GWSpec(kRs_vals, params1);
    
    // // add .plot() function for GWSpec (also include vw and alpha used in title)!!
    // plt::figure_size(800, 600);
    // plt::loglog(OmegaGW.k(), OmegaGW.P(), "k-");
    // plt::suptitle("vw = " + to_string_with_precision(vw) + ", alN = " + to_string_with_precision(alN));
    // plt::xlabel("K=kRs");
    // plt::ylabel("Omega_GW(K)");
    // plt::xlim(1e-3, 1e+3);
    // // plt::ylim(1e-5, 1e+0);
    // plt::grid(true);
    // plt::save("../GW_spectrum.png");
    

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    ProfilerStop();
    /******************************************************************/
    return 0;
}
