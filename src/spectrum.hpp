// spectrum.hpp
#pragma once

#include "PhaseTransition.hpp"

namespace Spectrum {

/**
 * @brief Calculates prefactor for GW power spectrum $\Omega_{GW}$
 *
 * @param csq The square of the speed of sound (in the broken phase)
 * @param T0 Temperature of the universe today
 * @param H0 Hubble constant today
 * @param g0 Number of dof today
 * @param gs Number of dof at beginning of FOPT
 *
 * @return Constant prefactor of eq (93) in arXiv:2308.12943
 */
double prefac(double csq, double T0, double H0, double g0, double gs);

/**
 * @brief Calculates prefactor for GW power spectrum $\Omega_{GW}$
 *
 * @param csq The square of the speed of sound (in the broken phase)
 * @param u Universe object (stores universe parameters at start of PT and present day)
 *
 * @return Constant prefactor of eq (93) in arXiv:2308.12943
 */
double prefac(double csq, const PhaseTransition::Universe &u);

} // namespace Spectrum