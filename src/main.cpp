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

    // test_FluidProfile();
    PhaseTransition::PTParams params;
    std::cout << "vpm=" << params.vpm()[0] << "," << params.vpm()[1] << std::endl;

    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;

    return 0;
}
