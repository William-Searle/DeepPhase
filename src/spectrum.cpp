// spectrum.cpp
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "maths_ops.hpp"
#include "PhaseTransition.hpp"
#include "hydrodynamics.hpp"
#include "spectrum.hpp"

/*
TO DO:
- update prefac to allow for non-bag model
- update prefac to do actual calculation of TGW, OmegaK_KK
- remove instances of std::pow when possible - it is slow
- change throw exception for P() and k() so that it uses Pvec() and kvec() when wrong one is called
- update Ekin to pass in Profile class (or maybe just PTParams?)
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

void PowerSpec::write(const std::string& filename) const {
    // if (typeid(data_).name ) {
    //     throw std::invalid_argument("Power spectrum is not a vector, cannot write to disk.")
    // }
    std::cout << "Writing power spectrum to disk... ";
    std::ofstream file(filename);
    file << "k,P\n";

    const auto k_vals = std::get<SpectrumVec>(data_).first;
    const auto P_vals = std::get<SpectrumVec>(data_).second;
    for (size_t i = 0; i < k_vals.size(); ++i) {
        file << k_vals[i] << "," << P_vals[i] << "\n";
    }
    file.close();
    std::cout << "Saved to " << filename << "!\n";

    return;
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

// tau_fin and tau_s could be passed in by params?
// double dlt(double tau_fin, double tau_s, double k, double p, double pt) {
//     auto integrand = [&](double tau_m, double k, double p, double pt) {
//         return ff(tau_m, p) * ff(tau_m, pt) * std::cos(k * tau_m)
//     };

    
// }

double ff(double tau_m, double k, double cs) {
    return std::cos(k * cs * tau_m); // for SSM -> NEED TO UPDATE THIS
}

double dtau_fin(double tau_fin, double tau_s) {
    return tau_fin - tau_s;
}

// get rid of std::pow
// PowerSpec Ekin(double k, double csq, double beta, double Rs, const std::string &nuc_type) {
//     // Integrand of Ekin
//     auto integrand = [&](double Ttilde) {
//         return Hydrodynamics::lifetime_dist(Ttilde, nuc_type) * power6(Ttilde) * Hydrodynamics::Ap_sq(Ttilde * k / beta, csq);
//     };

//     // Integration routine
//     boost::math::quadrature::gauss_kronrod<double, 15> integrator;
//     double f = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
//     f *= std::pow(k / M_PI, 2) / (2.0 * power6(beta) * std::pow(Rs, 3));
    
//     return PowerSpec(k, f);
// }

// PowerSpec Ekin(double k, Hydrodynamics::FluidProfile& prof) {
//     Ekin(k, prof.params().csq(), prof.params().beta(), prof.params().Rs(), prof.params().nuc_type())
// }

PowerSpec Ekin(double k, Hydrodynamics::FluidProfile& prof) {
    const auto csq = prof.params().csq();
    const auto beta = prof.params().beta();
    const auto Rs = prof.params().Rs();
    const auto nuc_type = prof.params().nuc_type();

    auto lt_dist = Hydrodynamics::lifetime_dist_func(nuc_type);

    // calculate and store Ap_sq vals:
    // std::vector<double> chi = logspace(0.01, 1000, 500);
    // const auto Apsq_vals = Hydrodynamics::Ap_sq();
    // Integrand of Ekin
    auto integrand = [&](double Ttilde) {
        return lt_dist(Ttilde) * power6(Ttilde) * Hydrodynamics::Ap_sq(Ttilde * k / beta, prof);
    };

    // Integration routine
    boost::math::quadrature::gauss_kronrod<double, 15> integrator;
    double f = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
    f *= std::pow(k / M_PI, 2) / (2.0 * power6(beta) * std::pow(Rs, 3));
    
    return PowerSpec(k, f);
}

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
    return prefac(csq, u.T0(), u.H0(), u.g0(), u.gs());
}

} // namespace Spectrum




// #include <iostream>
// #include <vector>
// #include <cmath>
// #include <functional>
// #include <boost/math/quadrature/gauss_kronrod.hpp>
// #include "CubicSpline.hpp" // Assuming your spline class is defined here

// constexpr double PI = 3.141592653589793;

// // Example input profile function (can be replaced with interpolated profiles)
// double example_v_ip(double xi) {
//     return std::exp(-xi * xi); // placeholder
// }

// double example_lambda_ip(double xi) {
//     return std::exp(-xi); // placeholder
// }

// // f(chi)
// double f_chi(double chi, const std::function<double(double)>& v_ip, double upper_limit = 20.0, double tol = 1e-6) {
//     auto integrand = [&](double xi) {
//         return xi * v_ip(xi) * std::sin(chi * xi);
//     };
//     double I = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, 0.0, upper_limit, tol);
//     return (4.0 * PI / chi) * I;
// }

// // Numerical derivative
// double numerical_derivative(const std::function<double(double)>& f, double chi, double delta = 1e-4) {
//     return (f(chi + delta) - f(chi - delta)) / (2.0 * delta);
// }

// int main() {
//     std::vector<double> chi_vals;
//     std::vector<double> eps_vals;

//     double chi_min = 1e-2;
//     double chi_max = 50.0;
//     double chi_step = 1.05; // logarithmic spacing
//     double cs2 = 1.0 / 3.0; // example value

//     for (double chi = chi_min; chi <= chi_max; chi *= chi_step) {
//         chi_vals.push_back(chi);

//         auto f = [&](double ch) { return f_chi(ch, example_v_ip); };
//         auto l = [&](double ch) { return f_chi(ch, example_lambda_ip); };

//         double f_p = numerical_derivative(f, chi);
//         double l_p = numerical_derivative(l, chi);

//         double eps = 0.25 * (f_p * f_p + cs2 * l_p * l_p);
//         eps_vals.push_back(eps);
//     }

//     CubicSpline<double> E_interp(chi_vals, eps_vals);

//     // Integrate over chi
//     auto E_integrand = [&](double chi) {
//         return E_interp(chi);
//     };

//     double final_result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(
//         E_integrand, chi_min, chi_max, 1e-6
//     );

//     std::cout << "Integrated E^(1)(chi): " << final_result << std::endl;
//     return 0;
// }
