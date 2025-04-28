// hydrodynamics.cpp

/*
TO DO:
- change short functions to inline functions and move to header
*/

#include <cmath>
#include <array>
#include <complex>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <functional>

#include "constants.hpp"

namespace Hydrodynamics {

double mu(double xi, double v) {
    return (xi - v) / (1.0 - xi * v);
}

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
std::vector<double> lifetime_dist2(const std::vector<double> &Ttilde, const std::string &nuc_type) {
    std::vector<double> dist;
    dist.reserve(Ttilde.size()); // allocates fixed memory to dist

    std::function<double(double)> lifetime_func;

    if (nuc_type == "exp") {
        lifetime_func = [](double Tt) { return std::exp(-Tt); };
    } else if (nuc_type == "sim"){
        lifetime_func = [](double Tt) {
            const auto exp_fac = -pow(Tt, 3) / 6.;
            return 0.5 * pow(Tt, 2) * std::exp(exp_fac);
        };
    }
    else {
        throw std::invalid_argument("Invalid nucleation type: " + nuc_type);
    }

    for (const auto &Tt : Ttilde)
        dist.push_back(lifetime_func(Tt));
    return dist;
}

// defines system of bubble profile DEs (not finished)
std::array<double,3> dfdv(const std::array<double,3>& xiw, double v, double csq) {
    if (xiw.size() != 3) {
        throw std::invalid_argument("xiw must be an array of the form {xi, w, T}");
    }
    const auto& [xi, w, T] = xiw;
    const auto dxidv  = (std::pow(mu(xi, v),2) / csq - 1) * (1 - v*xi) 
                        * (xi/2) * (1/v) / (1-std::pow(v,2));             // velocity profile
    const auto dwdv = (1 + 1/csq) * mu(xi, v) * w / (1 - std::pow(v,2)); // enthalpy profile
    const auto dTdv = mu(xi,v) * T / (1 - std::pow(v,2));                 // temperature profile

    // change this and xiw to vec?
    // set up template function so output matches data type of xiw?
    const std::array<double,3> derivs = {dxidv, dwdv, dTdv};

    return derivs;
}

// use eigen to solve for velocity and enthalpy profiles
// double v_w_prof(double xi) {
//     return {v_prof, w_prof};
// }

double v_prof(double chi) {
    // const auto vw = v_w_prof(chi);
    // return vw[0];
    return 1.0;
}

double w_prof(double chi) {
    // const auto vw = v_w_prof(chi);
    // return vw[1];
    return 1.0;
}

double la_prof(double chi) { // lambda = (rho-rhobar)/wbar
    return 1.0;
}

// not sure if needed
double prof_int_f(double chi) { // calculate integral in eq 30
    auto integrand = [&](double xi) {
        return v_prof(xi) * std::sin(chi * xi);
    };

    boost::math::quadrature::gauss_kronrod<double, 15> integrator;
    double f = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
    f *= 4.0 * PI / chi;

    return f;
}

double prof_int_f_der(double chi) {
    // could break up integrand into two integrals - effects efficiency?
    auto integrand = [&](double xi) {
        return v_prof(xi) * (xi * std::cos(chi * xi) - std::sin(chi * xi) / chi);
    };

    boost::math::quadrature::gauss_kronrod<double, 15> integrator;
    double f_der = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
    f_der *= 4.0 * PI / chi;

    return f_der;
}

double prof_int_l(double chi) {
    auto integrand = [&](double xi) {
        return xi * la_prof(xi) * std::sin(chi * xi);
    };

    boost::math::quadrature::gauss_kronrod<double, 15> integrator;
    double l = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
    l *= 4.0 * PI / chi;

    return l;
}

// not sure if needed
std::complex<double> Apm(std::string pm, double chi, const double cs) {
    double sgn;
    if (pm == "+") {
        sgn = 1;
    } else if (pm == "-") {
        sgn = -1;
    } else {
        throw std::invalid_argument("Invalid argument pm=" + pm + ". Requires '+' or '-'.");
    }

    std::complex<double> i(0.0, 1.0); // define elsewhere?
    const std::complex<double> AA = (-i/2.0) * (prof_int_f_der(chi) + sgn * i * cs * prof_int_l(chi));

    return AA;
}

// |A_+|^2
double Ap_sq(double chi, const double csq) {
    // pow(prof_int_f_der,2) calls the function twice, so this is more efficient
    const auto prof_int_1 = prof_int_f_der(chi);
    const auto prof_int_2 = prof_int_l(chi);
    return 0.25 * (prof_int_1 * prof_int_1) - csq * (prof_int_2 * prof_int_2);
}


} // namespace Hydrodynamics