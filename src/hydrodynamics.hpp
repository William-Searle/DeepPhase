// hydrodynamics.hpp
#pragma once

/*
TO DO:
- define a class NucModel that stores nuc_type (exp or sim), bubble lifetime distribution function v(Ttilde), deflag/det/hybrid type
*/

#include <array>
#include <complex>

#include "profile.hpp"

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

double lifetime_dist(double Ttilde, const std::string &nuc_type);

/**
 * @brief Calculates the bubble lifetime distribution for exponential/simultaneous nucleation
 *
 * @param Ttilde Ttilde = T \beta (normalised bubble lifetime vector)
 * @param nuc_type Nucleation type (exponential 'exp' or simultaneous 'sim')
 *
 * @return Bubble lifetime distribution
 */
std::vector<double> lifetime_dist2(const std::vector<double> &Ttilde, const std::string &nuc_type);

/**
 * @brief Calculates the integrated profile function 'f' (eq 30 in Pol, Procacci, Caprini (2024))
 *
 * @param chi ..
 *
 * @return f
 */
double prof_int_f(double chi, FluidProfile& prof);

double prof_int_f_der(double chi, FluidProfile& prof);

double prof_int_l(double chi, FluidProfile& prof);

std::complex<double> Apm(std::string pm, double chi, FluidProfile& prof);

double Ap_sq(double chi, FluidProfile& prof);

} // namespace Hydrodynamics