// spectrum.hpp
#pragma once

#include <vector>
#include <string>
#include <variant>
#include <algorithm>

#include "PhaseTransition.hpp"

/*
TO DO:
- use template in PowerSpec ctor -> only need one ctor
- destructor for PowerSpec to save memory somewhere?
- update operator overload to use type traits (for int, unsigned int, etc. scalar)
- PowerSpec class documentation
*/

namespace Spectrum {

class PowerSpec {
  public:
    using Spectrum = std::pair<double, double>;
    using SpectrumVec = std::pair<std::vector<double>, std::vector<double>>;

    // ctors - pass in momentum (k) and spectrum (P) since P is calculated differently for kinetic/GW
    PowerSpec(double k, double P);
    PowerSpec(const std::vector<double>& kvec, std::vector<double>& Pvec);

    double k() const; // Momentum
    const std::vector<double>& kvec() const;

    double P() const; // Power spectrum
    const std::vector<double>& Pvec() const;

    double max() const; // Max value of power spectrum
    bool is_scalar() const; // Move to private after testing?

    // Scalar arithmetic
    // Note: += changes an object in the class, but + creates a new object, 
    // so + is not a member function and we therefore define it as a friend 
    // (non-member function that can access private/protected parts of the object)
    friend PowerSpec operator+(const PowerSpec &spec, double scalar);
    friend PowerSpec operator+(double scalar, const PowerSpec &spec);
    PowerSpec &operator+=(double scalar);

    friend PowerSpec operator-(const PowerSpec &spec, double scalar);
    friend PowerSpec operator-(double scalar, const PowerSpec &spec);
    PowerSpec &operator-=(double scalar);

    friend PowerSpec operator*(const PowerSpec &spec, double scalar);
    friend PowerSpec operator*(double scalar, const PowerSpec &spec);
    PowerSpec &operator*=(double scalar);

    friend PowerSpec operator/(const PowerSpec& spec, double scalar);
    PowerSpec &operator/=(double scalar);
    // define "/" for PowerSpec on demoninator? Probably not needed...    

    // PowerSpec arithmetic
    friend PowerSpec operator+(const PowerSpec &spec1, const PowerSpec &spec2);
    PowerSpec &operator+=(PowerSpec& spec);

    friend PowerSpec operator-(const PowerSpec &spec1, const PowerSpec &spec2);
    PowerSpec &operator-=(PowerSpec &spec);

    friend PowerSpec operator*(const PowerSpec &spec1, const PowerSpec &spec2);
    PowerSpec &operator*=(PowerSpec &spec);

  private:
    std::variant<Spectrum, SpectrumVec> data_;

};
#include "spectrum.tpp"

  /**
   * @brief Calculates kinetic (velocity) power spectrum
   *
   * @param k Momentum
   * @param beta Inverse PT duration
   * @param Rs Characteristic length scale R_*
   * @param nuc_type Bubble nucleation method (exp or sim)
   *
   * @return Kinetic power spectrum
   */
  PowerSpec Ekin(double k, double csq, double beta, double Rs, const std::string &nuc_type);

  /**
   * @brief Calculates kinetic (velocity) power spectrum
   *
   * @param k Momentum
   * @param params Phase transition parameters
   *
   * @return Kinetic power spectrum
   */
  PowerSpec Ekin(double k, const PhaseTransition::PTParams &params);

  /**
   * @brief Calculates normalised kinetic power spectrum from Ekin(k)
   *
   * @param Ekin kinetic power spectrum
   *
   * @return Normalised kinetic power spectrum
   */
  PowerSpec zetaKin(PowerSpec Ekin);

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
 * @param u Universe parameters at start of PT and present day
 *
 * @return Constant prefactor of eq (93) in arXiv:2308.12943
 */
double prefac(double csq, const PhaseTransition::Universe &u);

} // namespace Spectrum