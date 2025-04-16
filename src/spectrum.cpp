// spectrum.cpp
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <iostream>

#include "maths_ops.hpp"
#include "constants.hpp"
#include "PhaseTransition.hpp"
#include "hydrodynamics.hpp"
#include "spectrum.hpp"

/*
TO DO:
- update prefac to allow for non-bag model
- update prefac to do actual calculation of TGW, OmegaK_KK
- remove instances of std::pow when possible - it is slow
- write linspace equivalent func
- write PowerSpec class
- add vec class to multiply/divide vectors by a scalar and take powers of it or some other way of doing this cleanly
*/

namespace Spectrum {

/***** PowerSpec class *****/

// Define ctors
PowerSpec::PowerSpec(double k, double P)
    : data_(Spectrum{k, P}) {}

PowerSpec::PowerSpec(const std::vector<double> &kvec, std::vector<double> &Pvec)
    : data_(SpectrumVec{kvec, Pvec}) {
        if (kvec.size() != Pvec.size()) {
            throw std::invalid_argument("PowerSpec: k and P vectors must be the same size!");
        }
    }

// Public functions
bool PowerSpec::is_scalar() const {
    return std::holds_alternative<Spectrum>(data_);
}

double PowerSpec::k() const {
    if (!is_scalar()) 
        throw std::runtime_error("PowerSpec: Called k() on vector PowerSpec, use kvec() instead!");
    return std::get<Spectrum>(data_).first;
}

double PowerSpec::P() const {
    if (!is_scalar())
        throw std::runtime_error("PowerSpec: Called P() on vector PowerSpec, use Pvec() instead!");
    return std::get<Spectrum>(data_).second;
}

const std::vector<double>& PowerSpec::kvec() const {
    if (is_scalar())
        throw std::runtime_error("PowerSpec: Called kvec() on scalar PowerSpec, use k() instead!");
    return std::get<SpectrumVec>(data_).first;
}

const std::vector<double>& PowerSpec::Pvec() const {
    if (is_scalar())
        throw std::runtime_error("PowerSpec: Called Pvec() on scalar PowerSpec, use P() instead!");
    return std::get<SpectrumVec>(data_).second;
}

double PowerSpec::max() const {
    if (is_scalar()) {
        std::cerr << "Warning: Called max() on scalar PowerSpec. Returning scalar value.\n";
        return P();
    }
    const auto &Pv = Pvec();
    return *std::max_element(Pv.begin(), Pv.end());
}

// PowerSpec [op] Scalar arithmetic
PowerSpec operator+(const PowerSpec& spec, double scalar) {
    return Spectrum::scalar_arith(spec, scalar, std::plus<>());
}

PowerSpec operator+(double scalar, const PowerSpec& spec) {
    return spec + scalar;
}

PowerSpec& PowerSpec::operator+=(double scalar) {
    *this = *this + scalar;
    return *this;
}

PowerSpec operator-(const PowerSpec& spec, double scalar) {
    return Spectrum::scalar_arith(spec, scalar, std::minus<>());
}

PowerSpec operator-(double scalar, const PowerSpec& spec) {
    return scalar + (-1.0) * spec;
}

PowerSpec& PowerSpec::operator-=(double scalar) {
    *this = *this - scalar;
    return *this;
}

PowerSpec operator*(const PowerSpec& spec, double scalar) {
    return Spectrum::scalar_arith(spec, scalar, std::multiplies<>());
}

PowerSpec operator*(double scalar, const PowerSpec& spec) {
    return spec * scalar;
}

PowerSpec& PowerSpec::operator*=(double scalar) {
    *this = *this * scalar;
    return *this;
}

PowerSpec operator/(const PowerSpec& spec, double scalar) {
    if (scalar == 0)
        throw std::invalid_argument("PowerSpec: Division by zero!");
    return spec * (1.0 / scalar);
}

PowerSpec& PowerSpec::operator/=(double scalar) {
    *this = *this / scalar;
    return *this;
}

// PowerSpec [op] PowerSpec arithmetic
PowerSpec operator+(const PowerSpec& spec1, const PowerSpec& spec2) {
    return Spectrum::spec_arith(spec1, spec2, std::plus<>());
}

PowerSpec& PowerSpec::operator+=(PowerSpec& spec) {
    *this = *this + spec;
    return *this;
}

PowerSpec operator-(const PowerSpec& spec1, const PowerSpec& spec2) {
    return Spectrum::spec_arith(spec1, spec2, std::minus<>());
}

PowerSpec& PowerSpec::operator-=(PowerSpec& spec) {
    *this = *this - spec;
    return *this;
}

PowerSpec operator*(const PowerSpec& spec1, const PowerSpec& spec2) {
    return Spectrum::spec_arith(spec1, spec2, std::multiplies<>());
}

PowerSpec& PowerSpec::operator*=(PowerSpec& spec) {
    *this = *this * spec;
    return *this;
}

/***************************/

// move to hydrodynamics or profile?
double lifetime_dist(double Ttilde, const std::string &nuc_type) {

    std::function<double(double)> lifetime_func;

    if (nuc_type == "exp") {
        return std::exp(-Ttilde);
    } else if (nuc_type == "sim") {
        const auto exp_fac = -pow(Ttilde, 3) / 6.;
        return 0.5 * pow(Ttilde, 2) * std::exp(exp_fac);
    } else {
        throw std::invalid_argument("Invalid nucleation type: " + nuc_type);
    }
}

// calculate as vector here rather than calling function at integration step since it is expensive
std::vector<double> lifetime_dist2(const std::vector<double> &Ttilde, const std::string &nuc_type)
{
    std::vector<double> dist;
    dist.reserve(Ttilde.size()); // allocates fixed memory to dist

    std::function<double(double)> lifetime_func;

    if (nuc_type == "exp")
    {
        lifetime_func = [](double Tt)
        {
            return std::exp(-Tt);
        };
    }
    else if (nuc_type == "sim")
    {
        lifetime_func = [](double Tt)
        {
            const auto exp_fac = -pow(Tt, 3) / 6.;
            return 0.5 * pow(Tt, 2) * std::exp(exp_fac);
        };
    }
    else
    {
        throw std::invalid_argument("Invalid nucleation type: " + nuc_type);
    }

    for (const auto &Tt : Ttilde)
        dist.push_back(lifetime_func(Tt));
    return dist;
}

// update when Ap_sq function finished
// PowerSpec Ekin(double k, double beta, double Rs, std::string &nuc_type) {    
//     auto integrand = [&](double Ttilde) {
//         // return lifetime_dist(Ttilde, nuc_type) * std::pow(Ttilde, 6) * Ap_sq(Ttilde * k / beta);
//         return lifetime_dist(Ttilde, nuc_type) * std::pow(Ttilde, 6);
//     };

//     boost::math::quadrature::tanh_sinh<double> integrator;
//     double f = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
//     f *= std::pow(k / PI, 2) / (2 * std::pow(beta, 6) * std::pow(Rs, 3));
    
//     return PowerSpec(k, f);
// }

// PowerSpec Ekin(double k, const PhaseTransition::PTParams &params) {
//     return Ekin(k, params.get_beta(), params.get_Rs(), params.get_nuc_type());
// }

PowerSpec zetaKin(PowerSpec Ekin) {
    return Ekin / Ekin.max();
}

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

} // namespace Spectrum