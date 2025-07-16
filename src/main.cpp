// main.cpp
#include <chrono>
#include <gperftools/profiler.h>

#include "tests.hpp"

int main() {
    /************************ CLOCK / PROFILER *************************/
    ProfilerStart("profile.out");
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    // Create and plot fluid profile
    // example_FluidProfile();

    // Create and plot kinetic power spectrum
    example_Kin_Spec();

    // Create and plot GW power spectrum
    example_GW_Spec();
    

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    ProfilerStop();
    /******************************************************************/
    return 0;
}