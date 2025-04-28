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

/*
TO DO:
- separate into different test files to test each function so I can compare to the python output
- use libraries?
- Input file where you can specify universe and pt constants
- test spectrum functions (Ekin, variants of delta, GW spectrum) behaviour in limits they discuss in paper
- define addition/subtraction/multiplication/division of vectors with vectors/scalars
*/

int main() {
    // const std::vector<double> Ttilde = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
    // const std::string nuc_type = "exp";
    // const auto distro = Spectrum::lifetime_dist(Ttilde, nuc_type);
    // for (int i = 1; i < Ttilde.size(); i++) {
    //     std::cout << Ttilde[i] << ", " << distro[i] << "\n";
    // }
    // const auto ti = std::chrono::high_resolution_clock::now();

    // const PhaseTransition::PTParams params;
    // const auto Ek = Spectrum::Ekin(0.5, params);

    // std::cout << "Ekin: (" << Ek.k() << "," << Ek.P() << ")" << std::endl;

    // const auto tf = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = tf - ti;
    // std::cout << "Timer: " << duration.count() << " s" << std::endl;
    const auto x = linspace(1e-10, 1e10, 100);
    for (double xx : x) {
        const auto Apsq = Hydrodynamics::Ap_sq(xx, 1./3.);
        std::cout << "(" << xx << "," << Apsq << ")" << "\n";
    }

    return 0;
}
