// hydrodynamics.cpp

/*
TO DO:
- change short functions to inline functions and move to header
*/

#include <cmath>
#include <array>
#include <complex>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <functional>

#include "constants.hpp"

namespace Hydrodynamics {

double mu(double xi, double v) {
    return (xi - v) / (1.0 - xi * v);
}

double getwow(double vp, double vm) {
    const auto wp = vp / (1. - std::pow(vp, 2));
    const auto wm = vm / (1. - std::pow(vm, 2));

    return wp/wm;
}

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

// not sure if needed
// std::complex<double> A(std::string pm, double chi, const double cs) {
//     double sgn;
//     if (pm == "+") {
//         sgn = 1;
//     } else if (pm == "-") {
//         sgn = -1;
//     } else {
//         throw std::invalid_argument("Invalid argument pm=" + pm + ". Requires '+' or '-'.");
//     }

//     std::complex<double> i(0.0, 1.0); // define elsewhere?
//     const std::complex<double> AA = (-i/2.0) * (prof_int_f_der(chi) + sgn * i * cs * prof_int_l(chi));

//     return AA;
// }

// std::vector<double> Ap_sq(std::vector<double> chi, const double csq) {
//     // pow(prof_int_f_der,2) calls the function twice, so this is more efficient
//     const auto prof_int_1 = prof_int_f_der(chi);
//     const auto prof_int_2 = prof_int_l(chi);
//     return 0.25 * (std::pow(prof_int_1,2) - csq * std::pow(prof_int_2,2));
// }

// double prof_int_f(double chi) {
//     // calculate integral in eq 30

//     auto integrand = [&](double xi) {
//         // return v(xi) * std::sin(chi * xi);
//         return std::sin(chi * xi);
//     };

//     boost::math::quadrature::tanh_sinh<double> integrator;
//     double f = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
//     f *= 4.0 * PI / chi;

//     return f;
// }

// double v_w_prof(double xi) {
//     // use eigen to solve for velocity and enthalpy profiles
//     return {v_prof, w_prof};
// }

// double prof_int_f_der(double chi) {
//     // calculate derivative of integral in eq 30
// }

// double prof_int_l(double chi) {
//     // calculate integral in eq 31

// }

} // namespace Hydrodynamics