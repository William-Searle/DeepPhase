// hydrodynamics.hpp
#pragma once

/*
TO DO:
- define a class NucModel that stores nuc_type (exp or sim), bubble lifetime distribution function v(Ttilde), deflag/det/hybrid type
*/

#include <array>

namespace Hydrodynamics {

/**
 * @brief Lorentz factor between wall frame and Universe frame
 *
 * @param xi Fluid velocity in wall frame
 * @param v Fluid velocity in universe frame
 *
 * @return Lorentz factor between wall frame and Universe frame
 */
double mu(double xi, double v);

/**
 * @brief The ratio of enthalpy of broken phase to the enthalpy of symmetric phase.
 *
 * @param vp Fluid velocity just infront of wall in wall frame
 * @param vm Fluid velocity just behind of wall in wall frame
 *
 * @return wp/wm
 */
double getwow(double vp, double vm);

/**
 * @brief The hydrodynamical equations for self-similar solution
 *
 * @param xi Fluid velocity in wall frame
 * @param v Fluid velocity in the universe frame
 * @param w Fluid enthalpy (more detail!)
 * @param T Fluid temperature (more detail!)
 * @param xiw Vector containing xi, w, T
 * @param csq The square of the speed of sound (in the broken phase)
 *
 * @return dxidv, dwdv, dTdv
 */
std::array<double,3> dfdv(const std::array<double,3>& xiw, double v, double cssq);

// double prof_int_f(double chi);

} // namespace Hydrodynamics