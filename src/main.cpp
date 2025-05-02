// main.cpp
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

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
    const auto ti = std::chrono::high_resolution_clock::now();

    // Hydrodynamics::state_type y0 = {0.1, 0.1, 0.1};
    // const Hydrodynamics::FluidProfile prof(y0, 1./3.);
    // prof.generate_streamplot_data();

    test_FluidProfile();

    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;

    return 0;
}
