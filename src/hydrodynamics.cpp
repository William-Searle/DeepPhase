// hydrodynamics.cpp

/*
TO DO:
- change short functions to inline functions and move to header
- check for domain errors when calling cublic spline in Ap_sq()
*/

#include <cmath>
#include <complex>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <functional>

#include "profile.hpp"

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

// not sure if needed
double prof_int_f(double chi, FluidProfile& prof) { // calculate integral in eq 30
    const auto v_prof = prof.v_prof();

    auto integrand = [&](double xi) {
        return v_prof(xi) * std::sin(chi * xi);
    };

    boost::math::quadrature::gauss_kronrod<double, 15> integrator;
    double f = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
    f *= 4.0 * M_PI / chi;

    return f;
}

double prof_int_f_der(double chi, FluidProfile& prof) {
    const auto v_prof = prof.v_prof();
    // could break up integrand into two integrals - effects efficiency?
    auto integrand = [&](double xi) {
        return v_prof(xi) * (xi * std::cos(chi * xi) - std::sin(chi * xi) / chi);
    };

    boost::math::quadrature::gauss_kronrod<double, 15> integrator;
    double f_der = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
    f_der *= 4.0 * M_PI / chi;

    return f_der;
}

double prof_int_l(double chi, FluidProfile& prof) {
    const auto la_prof = prof.la_prof();

    auto integrand = [&](double xi) {
        return xi * la_prof(xi) * std::sin(chi * xi);
    };

    boost::math::quadrature::gauss_kronrod<double, 15> integrator;
    double l = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
    l *= 4.0 * M_PI / chi;

    return l;
}

// not sure if needed
std::complex<double> Apm(std::string pm, double chi, FluidProfile& prof) {
    double sgn;
    if (pm == "+") {
        sgn = 1;
    } else if (pm == "-") {
        sgn = -1;
    } else {
        throw std::invalid_argument("Invalid argument pm=" + pm + ". Requires '+' or '-'.");
    }

    std::complex<double> i(0.0, 1.0); // define elsewhere?
    const std::complex<double> AA = (-i/2.0) * (prof_int_f_der(chi, prof) + sgn * i * std::sqrt(prof.csq()) * prof_int_l(chi, prof));

    return AA;
}

// |A_+|^2
// WARNING: need to be careful calling f and l integrals -> can cause domain errors since cubic spline has fixed domain (fix this here)
double Ap_sq(double chi, FluidProfile& prof) {
    // pow(prof_int_f_der(chi),2) calls the function twice, so this is more efficient
    const auto prof_int_1 = prof_int_f_der(chi, prof);
    const auto prof_int_2 = prof_int_l(chi, prof);
    return 0.25 * (prof_int_1 * prof_int_1) - prof.csq() * (prof_int_2 * prof_int_2);
}


} // namespace Hydrodynamics