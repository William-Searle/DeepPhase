// main.cpp
#include "hydrodynamics.hpp"
#include "profile.hpp"
#include "spectrum.hpp"
#include <iostream>
#include <iomanip>

int main() {
    // Initial values
    double v0 = 0.01;
    double w0 = 1.0;
    double v_end = 0.99;
    double eta = 1.0;
    double H = 1.0e-5;  // Hubble parameter (arbitrary units)
    double R = 1.0;     // Bubble radius
    double Uf = 1.0;    // Fluid velocity scale

    // Generate profile using RK4
    auto profile = generate_profile(v0, w0, v_end);

    // Compute kinetic energy
    double KE = computeKineticEnergy(profile.v_vals, profile.w_vals, v_end, eta);
    std::cout << "Kinetic energy: " << std::setprecision(6) << KE << "\n";

    // Generate transfer spectrum
    std::vector<double> k_vals;
    std::vector<double> D_vals;
    for (int i = 1; i <= 100; ++i) {
        double k = 0.1 * i;
        k_vals.push_back(k);
        D_vals.push_back(transferFunction(k * R));
    }

    // Compute GW spectrum
    auto omega = computeOmegaGW(k_vals, D_vals, H, R, Uf);

    // Print a few results
    std::cout << "\nSample GW Spectrum (k, Omega):\n";
    for (size_t i = 0; i < omega.size(); i += 10) {
        std::cout << "k = " << k_vals[i] << ", Omega = " << omega[i] << "\n";
    }

    return 0;
}
