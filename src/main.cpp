// main.cpp
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

#include "tests.hpp"
#include "constants.hpp"
#include "PhaseTransition.hpp"
#include "hydrodynamics.hpp"
#include "spectrum.hpp"
#include "maths_ops.hpp"
#include "profile.hpp"

/*
TO DO:
- separate into different test files to test each function so I can compare to the python output
- use libraries?
- Input file where you can specify universe and pt constants
- test spectrum functions (Ekin, variants of delta, GW spectrum) behaviour in limits they discuss in paper
- define addition/subtraction/multiplication/division of vectors with vectors/scalars
*/

int main() {
    const auto ti = std::chrono::high_resolution_clock::now();

    
    // const PhaseTransition::PTParams params;
    // const double k = 0.1;
    // const auto Ek = Spectrum::Ekin(k, params);

    // std::cout << Ek.k() << "," << Ek.P() << std::endl;

    // auto f = [](double x) { return Hydrodynamics::Ap_sq(x, 1./3.); };
    // plot_func(f, 0.1, 1., 1000);

    const std::vector<double> y0 = {0.1, 0.1, 0.1};
    FluidProfile::solve_prof(y0, 1./3.);
    FluidProfile::plot_velocity_profile();


    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;

    return 0;
}
