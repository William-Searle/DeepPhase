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
#include <omp.h>

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
        // std::cerr << "Warning: Called max() on scalar PowerSpec. Returning scalar value.\n";
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

/*** GW power spectrum ***/
// add option for inputing pRs_vals, Ttilde_vals and z_vals?
PowerSpec GWSpec(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params) {
    const auto Rs = params.Rs();
    const auto Rs_inv = 1.0 / Rs;

    const auto z_vals = linspace(-1.0, 1.0, 200);
    const auto nz = z_vals.size();

    const auto pRs_vals = linspace(1e-2, 1e+3, 200); // P = p*Rs
    const auto np = pRs_vals.size();

    std::vector<double> pRs2_vals; // keep here otherwise have to calculate for each k
    std::vector<double> p_vals; 
    for (const auto pRs : pRs_vals) {
        pRs2_vals.push_back(pRs * pRs);
        p_vals.push_back(pRs * Rs_inv);
    }

    // can probably get rid of this by modifying dlt to take in kRs and pRs rather than k and p
    std::vector<double> k_vals;
    for (const auto kRs : kRs_vals) {
        k_vals.push_back(kRs * Rs_inv);
    }

    // precompute normalised kinetic spectrum
    // create one zetaKin class object and define interpolating function to call inside loop
    // this avoids creating new PowerSpec objects in loop which is slow as hell
    const Hydrodynamics::FluidProfile profile(params);
    const auto zetaKin1 = zetaKin(pRs_vals, profile);
    const auto zetaKin1_P_vals = zetaKin1.Pvec();

    // precompute dlt
    const auto delta = dlt(k_vals, p_vals, z_vals, params);

    const auto nk = kRs_vals.size();
    std::vector<double> GW_P_vals(nk);

    #pragma omp parallel for
    for (int m = 0; m < nk; m++ ) {
        const auto kRs = kRs_vals[m];
        const auto k = kRs * Rs_inv;

        std::vector<std::vector<double>> integrand(np, std::vector<double>(nz));
        for (int i = 0; i < np; i++) {
            const auto p = p_vals[i];
            const auto pRs = pRs_vals[i];
            const auto pRs2 = pRs2_vals[i];

            const auto zetaKin1_fac = kRs * zetaKin1_P_vals[i] * pRs2; // kRs * zetaKin(pRs) * pRs^2

            for (int j = 0; j < nz; j++) {
                const auto z = z_vals[j];
                const auto pt = ptilde(k, p, z);
                const auto ptRs = pt * Rs;

                // calc 2nd kinetic spectrum (can't store earlier since pt=pt(k,p,z))
                // could calc ALL ptRs and then create only one Ek vector -> saves many ctor calls
                // this MUST be precomputed otherwise zetaKin is ALWAYS 1!!!
                const auto Ek = Ekin(ptRs, profile);
                const auto zetaKin2 = zetaKin(Ek);

                const auto dlta = delta[m][i][j];

                const auto ptRs4_inv = 1.0 / (ptRs * ptRs * ptRs * ptRs);
                const auto z_fac = 1.0 - z;
                const auto z_fac2 = z_fac * z_fac;

                integrand[i][j] = z_fac2 * ptRs4_inv * zetaKin1_fac * zetaKin2.P() * dlt;
                // integrand[i][j] = z_fac2 * ptRs4_inv * zetaKin1_fac * dlta;
            }
        }

        GW_P_vals[m] = simpson_2d_integrate(pRs_vals, z_vals, integrand);
    }

    return PowerSpec(kRs_vals, GW_P_vals);
}
/***************************/

kRs /int_0^{\infty} dP \int_{-1}^1 dz (1-z^2)^2 (P^2)/(Pt^4) zetaKin(P) zetaKin(Pt) dlt(k,p,z)


/*** dlt spectrum ***/
double ptilde(double k, double p, double z) {
    return std::sqrt(k*k - 2.0 * k * p * z + p*p);
}

double ff(double tau_m, double k, double cs) {
    return std::cos(k * cs * tau_m); // for SSM -> NEED TO UPDATE THIS
}

double dtau_fin(double tau_fin, double tau_s) {
    return tau_fin - tau_s;
}

// tau_vals is the same for each call of dlt() -> redundance when calling dlt() in a loop since it recalculates tau_m each time
// best to store dlt as a nk x np x npt tensor to avoid this
// WARNING: dlt takes in k, NOT K=kRs
std::vector<std::vector<std::vector<double>>> dlt(const std::vector<double>& k_vals, const std::vector<double>& p_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params) {
    const auto cs = std::sqrt(params.csq());

    const auto tau_s = params.tau_s();
    const auto tau_fin = params.tau_fin();

    const auto tau_vals = linspace(tau_s, tau_fin, 50);
    const auto n = tau_vals.size();

    // store tau_m = tau2 - tau1 and tau_sq_inv = 1/(tau1 * tau2) values
    // avoids repeated calculation in loops over k, p, z (much quicker!)
    std::vector<std::vector<double>> tau_m(n, std::vector<double>(n));
    std::vector<std::vector<double>> tau_sq_inv(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) { // tau1
        for (int j = 0; j < n; j++) { // tau2
            const auto tau1 = tau_vals[i];
            const auto tau2 = tau_vals[j];
            tau_m[i][j] = tau2 - tau1;
            tau_sq_inv[i][j] = 1.0 / (tau1 * tau2);
        }
    }

    // reserve memory for output tensor
    const auto nk = k_vals.size();
    const auto np = p_vals.size();
    const auto nz = z_vals.size();
    std::vector<std::vector<std::vector<double>>> result(nk, std::vector<std::vector<double>>(np, std::vector<double>(nz)));
    auto result2 = result;
    // bottleneck is subtraction in ptilde
    // better way to parallelise loops?
    // slightly slower initialising integrand inside z loop
    //  - could initialise outside of loop (with result) and give integrand a thread_id then sum all threads to make more efficient
    #pragma omp parallel for
    for (int kk = 0; kk < nk; kk++) {
        const auto k = k_vals[kk];
        for (int pp = 0; pp < np; pp++) {
            const auto p = p_vals[pp];
            for (int zz = 0; zz < nz; zz++) {
                const auto z = z_vals[zz];
                const auto pt = ptilde(k, p, z);
                std::vector<std::vector<double>> integrand(n, std::vector<double>(n));
                // integration routine
                for (int i = 0; i < n; i++) { // tau1
                    for (int j = 0; j < n; j++) { // tau2
                        const auto tau_minus = tau_m[i][j];
                        integrand[i][j] = ff(tau_minus, p, cs) * ff(tau_minus, pt, cs) * std::cos(k * tau_minus) * tau_sq_inv[i][j];
                    }
                }
                result[kk][pp][zz] = simpson_2d_integrate(tau_vals, tau_vals, integrand);
            }
        }
    }

    return result;
}

std::vector<std::vector<std::vector<double>>> dlt2(const std::vector<double>& k_vals, const std::vector<double>& p_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params) {
    const auto cs = std::sqrt(params.csq());

    const auto tau_s = params.tau_s();
    const auto tau_fin = params.tau_fin();

    const auto tau_vals = linspace(tau_s, tau_fin, 50);
    const auto n = tau_vals.size();

    // store tau_m = tau2 - tau1 and tau_sq_inv = 1/(tau1 * tau2) values
    // avoids repeated calculation in loops over k, p, z (much quicker!)
    std::vector<std::vector<double>> tau_m(n, std::vector<double>(n));
    std::vector<std::vector<double>> tau_sq_inv(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) { // tau1
        for (int j = 0; j < n; j++) { // tau2
            const auto tau1 = tau_vals[i];
            const auto tau2 = tau_vals[j];
            tau_m[i][j] = tau2 - tau1;
            tau_sq_inv[i][j] = 1.0 / (tau1 * tau2);
        }
    }

    // reserve memory for output tensor
    const auto nk = k_vals.size();
    const auto np = p_vals.size();
    const auto nz = z_vals.size();
    std::vector<std::vector<std::vector<double>>> result(nk, std::vector<std::vector<double>>(np, std::vector<double>(nz)));
    auto result2 = result;

    // splitting up ptilde and ff calcs like this doesn't really change efficiency
    #pragma omp parallel for
    for (int kk = 0; kk < nk; kk++) {
        const auto k = k_vals[kk];
        const auto k2 = k * k;
        for (int pp = 0; pp < np; pp++) {
            const auto p = p_vals[pp];
            const auto p2 = p * p;
            const auto kp = k * p;
            const auto kp2_sum = k2 + p2;
            const auto pcs = p * cs;
            for (int zz = 0; zz < nz; zz++) {
                const auto z = z_vals[zz];
                const auto pt = std::sqrt(kp2_sum - 2.0 * kp * z); // ptilde = sqrt(k^2 - 2kpz + p^2)
                const auto ptcs = pt * cs;

                std::vector<std::vector<double>> integrand(n, std::vector<double>(n));
                // integration routine
                for (int i = 0; i < n; i++) { // tau1
                    for (int j = 0; j < n; j++) { // tau2
                        const auto tau_minus = tau_m[i][j];
                        const auto ff_p = std::cos(pcs * tau_minus); // for SSM only
                        const auto ff_pt = std::cos(ptcs * tau_minus); // for SSM only

                        integrand[i][j] = ff_p * ff_pt * std::cos(k * tau_minus) * tau_sq_inv[i][j];
                    }
                }

                result[kk][pp][zz] = simpson_2d_integrate(tau_vals, tau_vals, integrand);
            }
        }
    }

    return result;
}
/***************************/

/*** Kinetic spectrum ***/
PowerSpec Ekin(double kRs, const Hydrodynamics::FluidProfile& prof) {
    const auto csq = prof.params().csq(); // better way of doing this?
    const auto beta = prof.params().beta();
    const auto Rs = prof.params().Rs();
    const auto nuc_type = prof.params().nuc_type();

    const auto k = kRs / Rs; // make this better so less * and / in function call
    const auto fac1 = beta / k;
    const auto fac2 = k / M_PI;
    const auto fac3 = fac2 * fac2 / (2.0 * power3(beta * beta * Rs));

    auto lt_dist = Hydrodynamics::lifetime_dist_func(nuc_type);

    // define Ttilde from chi = Ttilde * k / beta (makes calling Apsq simpler)
    // Ap_sq = inf at 0
    const auto chi_vals = linspace(0.1, 40.0, 200); // bad to hard code?
    const auto n = chi_vals.size();

    const auto Apsq = Hydrodynamics::Ap_sq(chi_vals, prof);

    // see if index-based loop below actually improves performance at all
    std::vector<double> Ttilde_vals(n), integrand(n);
    for (int i = 0; i < n; i++) {
        const auto Ttilde = fac1 * chi_vals[i];
        Ttilde_vals[i] = Ttilde;
        integrand[i] = fac3 * lt_dist(Ttilde) * power6(Ttilde) * Apsq[i];
    }

    double P = simpson_integrate(Ttilde_vals, integrand);

    return PowerSpec(k, P);
}

// takes in K=k*Rs
PowerSpec Ekin(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& prof) {
    const auto csq = prof.params().csq();
    const auto beta = prof.params().beta();
    const auto Rs = prof.params().Rs();
    const auto nuc_type = prof.params().nuc_type();

    auto lt_dist = Hydrodynamics::lifetime_dist_func(nuc_type);

    // define Ttilde from chi = Ttilde * k / beta (makes calling Apsq simpler)
    // using K = k * Rs below
    // Ap_sq = inf at 0
    const auto chi_vals = logspace(0.1, 40.0, 200); // bad to hard code?
    const auto n = chi_vals.size();

    const auto Apsq = Hydrodynamics::Ap_sq(chi_vals, prof);
    std::vector<double> P_vals;

    const auto fac1 = beta * Rs * Rs / (2.0 * M_PI * M_PI);
    // #pragma omp parallel for
    for (const auto kRs : kRs_vals) {
        const auto kRs_inv = 1.0 / kRs;
        const auto fac2 = fac1 * power(kRs_inv, 5);
        const auto fac3 = beta * Rs * kRs_inv;

        std::vector<double> integrand(n);
        for (int i = 0; i < n; i++) {
            const auto chi = chi_vals[i];
            integrand[i] = fac2 * lt_dist(fac3 * chi) * power(chi, 6) * Apsq[i];
            
        }

        const auto P = simpson_integrate(chi_vals, integrand);
        P_vals.push_back(P);
    }

    return PowerSpec(kRs_vals, P_vals);
}

PowerSpec zetaKin(const PowerSpec& Ekin) {
    const auto Ekin_max = Ekin.max();
    if (Ekin_max == 0.0) {
        throw std::runtime_error("Division by zero in zetaKin from Ekin.max() = 0");
    } else if (isnan(Ekin_max)) {
        throw std::runtime_error("In zetaKin: Ekin.max() = nan");
    }
    return Ekin / Ekin_max;
}

PowerSpec zetaKin(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& prof) {
    const auto Ek = Ekin(kRs_vals, prof);
    return zetaKin(Ek);
}
/***************************/

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