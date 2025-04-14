// spectrum.cpp
#include <cmath>
#include "PhaseTransition.hpp"

/*
TO DO:
- update prefac to allow for non-bag model
- update prefac to do actual calculation of TGW, OmegaK_KK
*/

namespace Spectrum {

double prefac(double csq, double T0, double H0, double g0, double gs) {
    // Gamma = w/e = 1+p/e (ratio of enthalpy to energy density)
    const auto Gamma = 1 + csq; // for bag model

    // these are vals used in paper - include actual calculation here
    const auto TGW = 1.0; // Transfer function (eq 13)
    const auto OmegaK_KK = 1e-4; // Omega_K / KK (eq 42)

    return 3 * std::pow(Gamma,2) * TGW * std::pow(OmegaK_KK,2);
}

double prefac(double csq, const PhaseTransition::Universe &u) {
    return prefac(csq, u.get_T0(), u.get_H0(), u.get_g0(), u.get_gs());
}

// std::vector<double> tildeDlt()

} // namespace Spectrum