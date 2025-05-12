// hydrodynamics.cpp

/*
TO DO:
- change short functions to inline functions and move to header
- check for domain errors when calling cublic spline in Ap_sq()
- could change f and l integrations to use profile.v_vals() rather than calling profile.v_prof(xi) every time (might be faster?)
    - i'm not sure if i would even need an interpolator if this is the case
- NOTE: on profile integrals f and l, xi is bound by 0<xi<1 so integration to infinity not strictly necessary (int from 0->1 and 0->inf agree to numerical precision)
    - might be safer to keep int to inf, but a bit quicker to just integrate to 1
- construct spline interpolation for f' and l integrands and use boost.math or integrate over a vector of integrand values using custom simpson/trapezoidal integrator?
    - for f' and l, they seem to be quite well behaved so both would be okay, but Ap_sq is highly oscillatory so i could reuse a custom integrator here (fitting spline to highly oscillatory function not good)
    - if adaptive integration (changing step-size for different parts of the func) is needed, calc vector of integrand values and construct a lambda from this somehow to use boost.math (avoids spline)
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
        const auto exp_fac = - Ttilde * Ttilde * Ttilde / 6.0;
        return 0.5 * Ttilde * Ttilde * std::exp(exp_fac);
    } else {
        throw std::invalid_argument("Invalid nucleation type: " + nuc_type);
    }
}

std::function<double(double)> lifetime_dist_func(const std::string &nuc_type) {
    // Return a lambda function that computes the lifetime distribution for the specified nucleation type
    if (nuc_type == "exp") {
        return [](double Ttilde) -> double {
            return std::exp(-Ttilde);
        };
    } else if (nuc_type == "sim") {
        return [](double Ttilde) -> double {
            const auto exp_fac = - Ttilde * Ttilde * Ttilde / 6.0;
            return 0.5 * Ttilde * Ttilde * std::exp(exp_fac);
        };
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

/* profile integrals (single value) */
double prof_int_f_der(double chi, FluidProfile& prof) {
    const auto v_prof = prof.v_prof();

    auto integrand = [&](double xi) {
        return v_prof(xi) * (xi * std::cos(chi * xi) - std::sin(chi * xi) / chi);
    };

    boost::math::quadrature::gauss_kronrod<double, 15> integrator;
    double f_der = integrator.integrate(integrand, 0.0, std::numeric_limits<double>::infinity());
    // double f_der = integrator.integrate(integrand, 0.0, 1000.0);
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
    // double l = integrator.integrate(integrand, 0.0, 10000.0);
    l *= 4.0 * M_PI / chi;

    return l;
}

/* profile integrals (vector) */
std::pair<std::vector<double>, std::vector<double>> prof_ints_fl(std::vector<double> chi_vals, FluidProfile& prof) {
    const auto xi_vals = prof.xi_vals();
    const auto v_vals = prof.v_vals();
    const auto la_vals = prof.la_vals();
    
    const auto n = xi_vals.size();
    const auto m = chi_vals.size();
    const auto fac = 4.0 * M_PI;

    // allocate memory for integrand/result
    std::vector<double> fd_integrand(n), fd(m); 
    std::vector<double> l_integrand(n), l(m);

    // integrand/integral evaluation
    for (int j = 0; j < m; j++) { // chi
        const auto chi = chi_vals[j];
        const auto inv_chi = 1.0 / chi;     

        for (int i = 0; i < n; i++) { // xi
            const auto xi = xi_vals[i];
            const auto v_prof = v_vals[i];
            const auto la_prof = la_vals[i];
            
            const auto chi_xi = chi * xi;
            const auto sin_cx = std::sin(chi_xi);
            const auto cos_cx = std::cos(chi_xi);

            fd_integrand[i] = fac * inv_chi * v_prof * (xi * cos_cx - sin_cx * inv_chi);
            l_integrand[i] = fac * inv_chi * xi * la_prof * sin_cx;
        }

        // much faster than interpolating function + boost integrator
        fd[j] = simpson_integrate(xi_vals, fd_integrand);
        l[j] = simpson_integrate(xi_vals, l_integrand);
    }
    
    return {fd, l}; // {f', l}
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
    const std::complex<double> AA = (-i/2.0) * (prof_int_f_der(chi, prof) + sgn * i * std::sqrt(prof.params().csq()) * prof_int_l(chi, prof));

    return AA;
}

// |A_+|^2
double Ap_sq(double chi, FluidProfile& prof) {
    const auto prof_int_1 = prof_int_f_der(chi, prof);
    const auto prof_int_2 = prof_int_l(chi, prof);

    return 0.25 * (prof_int_1 * prof_int_1) + prof.params().csq() * (prof_int_2 * prof_int_2);
}

// precomputes f' and l for all chi for efficiently
std::vector<double> Ap_sq(std::vector<double> chi_vals, FluidProfile& prof) {
    const auto csq = prof.params().csq();
    const auto [fd_int, l_int] = prof_ints_fl(chi_vals, prof);
    const auto m = chi_vals.size();

    std::vector<double> Apsq(m);
    for (int j = 0; j < m; j++) {
        const auto f = fd_int[j];
        const auto l = l_int[j];
        Apsq[j] = 0.25 * (f*f + csq * l*l);
    }

    return Apsq;
}

} // namespace Hydrodynamics