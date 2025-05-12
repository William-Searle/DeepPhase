// hydrodynamics.hpp
#pragma once

/*
TO DO:
- 
*/

#include <array>
#include <complex>

#include "profile.hpp"

namespace Hydrodynamics {

/**
 * @brief Computes the Lorentz factor between the wall frame and the universe frame.
 *
 * @param xi Fluid velocity in the wall frame.
 * @param v Fluid velocity in the universe frame.
 * 
 * @return Lorentz factor (mu).
 */
double mu(double xi, double v);

/**
 * @brief Computes the bubble lifetime distribution at the normalized time Ttilde.
 *
 * @param Ttilde Normalized time (T * β).
 * @param nuc_type Nucleation type: "exp" for exponential or "sim" for simultaneous.
 *
 * @return Lifetime distribution.
 */
double lifetime_dist(double Ttilde, const std::string &nuc_type);

/**
 * @brief Computes the bubble lifetime distribution over a vector of normalized times Ttilde.
 *
 * @param Ttilde Normalized time (T * β).
 * @param nuc_type Nucleation type: "exp" for exponential or "sim" for simultaneous.
 * 
 * @return Lifetime distribution.
 */
std::vector<double> lifetime_dist2(const std::vector<double> &Ttilde, const std::string &nuc_type);

std::function<double(double)> lifetime_dist_func(const std::string &nuc_type);

/**
 * @brief Computes the integrated profile function f(χ) defined in Eq. (30) of Pol, Procacci, Caprini (2024).
 *
 * @param chi χ=k*T_n (k = momentum, T_n = lifetime of n'th bubble)
 * @param prof FluidProfile object used for interpolation and evaluation.
 * 
 * @return Value of the integrated profile function f.
 */
double prof_int_f(double chi, FluidProfile& prof);

/**
 * @brief Computes the first derivative of the integrated profile function f(χ) defined in Eq. (30) of Pol, Procacci, Caprini (2024).
 *
 * @param chi χ=k*T_n (k = momentum, T_n = lifetime of n'th bubble)
 * @param prof FluidProfile object used for interpolation and evaluation.
 * 
 * @return Value of the integrated profile function f.
 */
double prof_int_f_der(double chi, FluidProfile& prof);

std::vector<double> prof_int_f_der(std::vector<double> chi_vals, FluidProfile& prof);

/**
 * @brief Computes the integrated profile function l(χ) defined in Eq. (31) of Pol, Procacci, Caprini (2024).
 *
 * @param chi χ=k*T_n (k = momentum, T_n = lifetime of n'th bubble)
 * @param prof FluidProfile object used for interpolation and evaluation.
 * 
 * @return Value of the integrated profile function l.
 */
double prof_int_l(double chi, FluidProfile& prof);

std::pair<std::vector<double>, std::vector<double>> prof_ints_fl(std::vector<double> chi_vals, FluidProfile& prof);

/**
 * @brief Computes the complex amplitude A₊ or A₋ defined in Eq. (29) of Pol, Procacci, Caprini (2024).
 *
 * @param pm String indicating plus or minus branch ("+" or "-").
 * @param chi χ=k*T_n (k = momentum, T_n = lifetime of n'th bubble)
 * @param prof FluidProfile object.
 * 
 * @return Complex amplitude A₊ or A₋.
 */
std::complex<double> Apm(std::string pm, double chi, FluidProfile& prof);

/**
 * @brief Computes |A₊(χ)|²
 *
 * @param chi χ=k*T_n (k = momentum, T_n = lifetime of n'th bubble)
 * @param prof Fluid profile object.
 * 
 * @return Squared modulus |A₊(χ)|².
 */
double Ap_sq(double chi, FluidProfile& prof);
std::vector<double> Ap_sq(std::vector<double> chi_vals, FluidProfile& prof);

} // namespace Hydrodynamics