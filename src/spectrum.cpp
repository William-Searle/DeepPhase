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
#include <chrono>
// #include <Eigen/Dense>

#include "maths_ops.hpp"
#include "PhaseTransition.hpp"
#include "hydrodynamics.hpp"
#include "spectrum.hpp"
#include "physics.hpp"

/*
TO DO:
- update prefac to allow for non-bag model
- update prefac to do actual calculation of TGW, OmegaK_KK
- remove instances of std::pow when possible - it is slow
- change throw exception for P() and k() so that it uses Pvec() and kvec() when wrong one is called
- update Ekin to pass in Profile class (or maybe just PTParams?)
- implement adaptive step-size in Ekin integration (and dlt later too)
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

CubicSpline<double> PowerSpec::interpolate() const {
    if (is_scalar()) {
        throw std::runtime_error("PowerSpec: Cannot interpolate scalar spectrum!");
    }
    return CubicSpline(kvec(), Pvec());
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

    const auto z_vals = linspace(-1.0, 1.0, 200); // logspace gives nan over this domain
    const auto nz = z_vals.size();

    const auto pRs_vals = logspace(1e-2, 1e+3, 200); // P = p*Rs
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

    /********** precompute normalised kinetic spectrum **********/
    // NOTE: this is currently a bit buggy and domain of interpolating function requires fine tuning depending on range of k,p,z
    // zetaKin(ptRs) can't be precomputed since ptRs = ptRs(k,p,z) -> use interpolator (much faster than constructing PowerSpec objects inside loops)

    // calc temp ptRs vals for interpolating func
    const auto pRs_max = pRs_vals.back();
    const auto kRs_max = kRs_vals.back();
    const auto ptRs_max = 10.0 * (kRs_max + pRs_max); // max of pt=sqrt(k^2-2kpz+p^2)
    const auto ptRs_min = 1e-7; // not sure how to choose best min val - update later

    const auto ptRs_vals_tmp = logspace(ptRs_min, ptRs_max, 2.0*np);

    const Hydrodynamics::FluidProfile profile(params); // generate fluid profile

    const auto zk_pRs_spec = zetaKin(pRs_vals, profile);
    const auto zk_pRs_vals = zk_pRs_spec.Pvec(); // store zetaKin(pRs) vals (quicker than calling interpolator)

    const auto zk_ptRs_spec = zetaKin(ptRs_vals_tmp, profile);    
    const auto zk_ptRs_interp = zk_ptRs_spec.interpolate(); // interpolating function for zetaKin(ptRs)
    /************************************************************/

    // precompute dlt
    const int nt = 50;
    const auto delta = dlt(nt, k_vals, p_vals, z_vals, params);

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

            const auto zk_pRs_fac = kRs * zk_pRs_vals[i] * pRs2; // kRs * zetaKin(pRs) * pRs^2

            for (int j = 0; j < nz; j++) {
                const auto z = z_vals[j];
                const auto pt = ptilde(k, p, z);
                const auto ptRs = pt * Rs;

                const auto ptRs4_inv = 1.0 / (ptRs * ptRs * ptRs * ptRs);
                const auto z_fac = 1.0 - z;
                const auto z_fac2 = z_fac * z_fac;

                // neccessary to store these? only called once so maybe not
                const auto zk_ptRs = zk_ptRs_interp(ptRs); // see ptRs_min/max if domain issues arise
                // const auto zk_ptRs = 1.0;
                const auto dlta = delta[m][i][j];

                integrand[i][j] = z_fac2 * ptRs4_inv * zk_pRs_fac * zk_ptRs * dlta;
            }
        }

        GW_P_vals[m] = simpson_2d_integrate(pRs_vals, z_vals, integrand);
    }

    return PowerSpec(kRs_vals, GW_P_vals);
}
/***************************/

/*** dlt spectrum ***/
double ptilde(double k, double p, double z) {
    const auto arg = k*k - 2.0 * k * p * z + p*p;
    return std::sqrt(std::max(arg, 0.0)); // gives 0 for large k,p (FIND BETTER FIX. THIS IS BAD IN CASE ARG <0 FROM BAD K,P)
}

double ff(double tau_m, double kcs) {
    // kcs = k*cs -> ff called this way to make dlt faster
    return std::cos(kcs * tau_m); // for SSM -> NEED TO UPDATE THIS
}

double dtau_fin(double tau_fin, double tau_s) {
    return tau_fin - tau_s;
}

// tau_vals is the same for each call of dlt() -> redundance when calling dlt() in a loop since it recalculates tau_m each time
// best to store dlt as a nk x np x npt tensor to avoid this
// WARNING: dlt takes in k, NOT K=kRs
std::vector<std::vector<std::vector<double>>> dlt(const int nt, const std::vector<double>& k_vals, const std::vector<double>& p_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params) {
    /***************************** CLOCK ******************************/
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    const auto cs = std::sqrt(params.csq());

    const auto tau_s = params.tau_s();
    const auto tau_fin = params.tau_fin();

    // integrand becomes very large for small tau -> use logspace for accuracy of integration
    // const auto nt = 50;
    const auto tau_vals = logspace(tau_s, tau_fin, nt);
    const auto ntsq = nt * nt;

    const auto nk = k_vals.size();
    const auto np = p_vals.size();
    const auto nz = z_vals.size();

    // store tau_m = tau2 - tau1 and tau_sq_inv = 1/(tau1 * tau2) values
    // avoids repeated calculation in loops over k, p, z (much quicker!)
    std::vector<double> tau_m(ntsq);
    std::vector<double> tau_sq_inv(ntsq);
    #pragma omp parallel for
    for (int idx = 0; idx < ntsq; idx++) {
        const int i = idx / nt; // idx = i * nt + j;
        const int j = idx % nt;

        const auto tau1 = tau_vals[i];
        const auto tau2 = tau_vals[j];

        tau_m[idx] = tau2 - tau1;
        tau_sq_inv[idx] = 1.0 / (tau1 * tau2);
    }

    // fill ff (reduces redundancy)
    std::vector<std::vector<double>> ff1_cache(ntsq, std::vector<double>(np));
    std::vector<std::vector<double>> ff3_cache(ntsq, std::vector<double>(nk));
    #pragma omp parallel for
    for (int idx = 0; idx < ntsq; idx++) {
        const auto tau_minus = tau_m[idx];
        // fill ff3
        for (int kk = 0; kk < nk; kk++) {
            const auto k = k_vals[kk];
            ff3_cache[idx][kk] = std::cos(k * tau_minus);
        }
        // fill ff1
        for (int pp = 0; pp < np; pp++) {
            const auto p = p_vals[pp];
            ff1_cache[idx][pp] = ff(tau_minus, p*cs); // redundancy computing p*cs here - change?
        }
    }

    // reserve memory for integration
    std::vector<std::vector<std::vector<double>>> result(nk, std::vector<std::vector<double>>(np, std::vector<double>(nz)));
    std::vector<std::vector<double>> integrands(omp_get_max_threads(), std::vector<double>(ntsq));

    // precompute weights for integration
    const auto weights = precompute_simpson_weights_2d(tau_vals, tau_vals);
    const auto Ax_weights = weights.Ax_weights;
    const auto Ay_weights = weights.Ay_weights;
    const auto dx = weights.dx;
    const auto dy = weights.dy;

    // integration routine
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<double>& integrand = integrands[tid]; // local thread copy of integrand

        #pragma omp for
        for (int idx = 0; idx < nk * np * nz; idx++) {
            const int kk = idx / (np * nz);
            const int pp = (idx / nz) % np;
            const int zz = idx % nz;
        // #pragma omp for collapse(3) schedule(dynamic)
        // for (int kk = 0; kk < nk; kk++)
        // for (int pp = 0; pp < np; pp++)
        // for (int zz = 0; zz < nz; zz++) {
            const auto k = k_vals[kk];
            const auto p = p_vals[pp];
            const auto z = z_vals[zz];

            if (std::isnan(z)) {
                throw std::runtime_error("bad z");
            }

            const auto pt = ptilde(k, p, z); // collapsing loops a lot quicker than breaking up ptilde calc so redundancy here is ok!
            const auto ptcs = pt * cs;

            // integration routine
            for (int i = 0; i < ntsq; i++) {
                const auto tau_minus = tau_m[i];
                const auto ff3 = ff3_cache[i][kk];
                const auto ff1 = ff1_cache[i][pp];
                const auto ff2 = ff(tau_minus, ptcs);

                integrand[i] = ff1 * ff2 * ff3 * tau_sq_inv[i];
            }

            // result[kk][pp][zz] = simpson_2d_nonuniform_flat(tau_vals, tau_vals, integrand);
            result[kk][pp][zz] = simpson_2d_nonuniform_flat_weighted(tau_vals, tau_vals, integrand, Ax_weights, Ay_weights, dx, dy);
        }
    }

    /***************************** CLOCK ******************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer (dlt): " << duration.count() << " s" << std::endl;
    /******************************************************************/

    return result;
}

// adaptive integration (very slow)
std::vector<std::vector<std::vector<double>>> dlt_adaptive(const std::vector<double>& k_vals, const std::vector<double>& p_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params) {
    /***************************** CLOCK ******************************/
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    const auto cs = std::sqrt(params.csq());

    const auto tau_s = params.tau_s();
    const auto tau_fin = params.tau_fin();

    const auto nk = k_vals.size();
    const auto np = p_vals.size();
    const auto nz = z_vals.size();

    // lambda function for integrand
    auto integrand = [cs](const double k, const double p, const double pt, double tau1, double tau2) {
        const auto tau_m = tau2 - tau1;

        const auto ff1 = std::cos(p * cs * tau_m);
        const auto ff2 = std::cos(pt * cs * tau_m);
        const auto ff3 = std::cos(k * tau_m);

        return ff1 * ff2 * ff3 / (tau1 * tau2);
    };

    // reserve memory for integration
    std::vector<std::vector<std::vector<double>>> result(nk, std::vector<std::vector<double>>(np, std::vector<double>(nz)));
    #pragma om parallel
    {
        #pragma omp for collapse(3) schedule(dynamic)
        for (int kk = 0; kk < nk; kk++)
        for (int pp = 0; pp < np; pp++)
        for (int zz = 0; zz < nz; zz++) {
            const auto k = k_vals[kk];
            const auto p = p_vals[pp];
            const auto z = z_vals[zz];

            const auto pt = ptilde(k, p, z); // collapsing loops a lot quicker than breaking up ptilde calc so redundancy here is ok!

            auto integrand_temp = [k, p, pt, &integrand](double tau1, double tau2) { // eval integrand at each k,p,z
                return integrand(k, p, pt, tau1, tau2);
            };

            result[kk][pp][zz] = adaptive_simpson_2d(integrand_temp, tau_s, tau_fin, tau_s, tau_fin, 1e-4, 5);
        }
    }

    /***************************** CLOCK ******************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer (dlt-adaptive): " << duration.count() << " s" << std::endl;
    /******************************************************************/

    return result;
}

// change call to Ci(), Si() to interpolating functions for these integrals (avoids calculating them on the fly and can just use function call)
std::vector<std::vector<std::vector<double>>> dlt_SSM(const std::vector<double>& k_vals, const std::vector<double>& p_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params) {
    /***************************** CLOCK ******************************/
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    const auto cs = std::sqrt(params.csq());

    const auto tau_s = params.tau_s();
    const auto tau_fin = params.tau_fin();

    const auto nk = k_vals.size();
    const auto np = p_vals.size();
    const auto nz = z_vals.size();

    // reserve memory for integration
    std::vector<std::vector<std::vector<double>>> result(nk, std::vector<std::vector<double>>(np, std::vector<double>(nz)));
    const std::vector<double> sum_vals = {-1.0, 1.0};

    #pragma omp parallel
    {
        const auto num_threads = omp_get_num_threads();
        const auto thread_id = omp_get_thread_num();

        // std::vector<double> dlt_temp(num_threads, 0.0);

        #pragma omp for collapse(3) schedule(dynamic)
        for (int kk = 0; kk < nk; kk++)
        for (int pp = 0; pp < np; pp++)
        for (int zz = 0; zz < nz; zz++) {
            const auto k = k_vals[kk];
            const auto p = p_vals[pp];
            const auto z = z_vals[zz];

            const auto pt = ptilde(k, p, z);
            auto dlt_temp = 0.0;

            for (const auto m : sum_vals) { // fill pmn, dlt
                const auto pmn_1 = (p + m * pt) * cs;
                for (const auto n : sum_vals) {
                    const auto pmn = pmn_1 + n * k;
                    const auto x1 = pmn * tau_fin;
                    const auto x2 = pmn * tau_s;

                    const auto dSi = Si(x1) - Si(x2);
                    const auto dCi = Ci(x1) - Ci(x2);

                    if (isnan(dSi)) {
                        std::cout << "dSi=nan for x1=" << x1 << ", x2=" << x2 << "\n";
                    }
                    // } else if (isnan(dCi)) {
                    //     std::cout << "dCi=nan for x1=" << x1 << ", x2=" << x2 << "\n";
                    // }

                    // dlt_temp[thread_id] += 0.25 * (dCi * dCi + dSi * dSi);
                    dlt_temp += 0.25 * (dCi * dCi + dSi * dSi);
                }
            }

            // result[kk][pp][zz] = dlt_temp[thread_id];
            result[kk][pp][zz] = dlt_temp;
        }
    }

    /***************************** CLOCK ******************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer (dlt_SSM): " << duration.count() << " s" << std::endl;
    /******************************************************************/

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
    const auto chi_vals = logspace(1e-2, 100, 100); // bad to hard code?
    const auto n = chi_vals.size();

    const auto Apsq = Hydrodynamics::Ap_sq(chi_vals, prof);

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

    const auto nk = kRs_vals.size();
    std::vector<double> P_vals(nk);
    // std::vector<double> P_vals;


    // define Ttilde from chi = Ttilde * k / beta (makes calling Apsq simpler)
    // using K = k * Rs below
    /*
    NOTE:
    - Ap_sq = inf at 0
    - chi_vals = logspace(1e-3, 3000, 5000) gives good convergence
    */
    const auto chi_vals = logspace(1e-3, 3000, 5000); // bad to hard code?
    const auto n = chi_vals.size();

    const auto Apsq = Hydrodynamics::Ap_sq(chi_vals, prof);

    const auto fac1 = beta * Rs * Rs / (2.0 * M_PI * M_PI);

    std::vector<std::vector<double>> integrands(omp_get_max_threads(), std::vector<double>(n));
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<double>& integrand = integrands[tid]; // local thread copy of integrand

        #pragma omp for
        for (int kk = 0; kk < nk; kk++) {
            const auto kRs = kRs_vals[kk];
            const auto kRs_inv = 1.0 / kRs;

            const auto fac2 = fac1 * power(kRs_inv, 5);
            const auto fac3 = beta * Rs * kRs_inv;

            for (int i = 0; i < n; i++) {
                const auto chi = chi_vals[i];
                integrand[i] = fac2 * lt_dist(fac3 * chi) * power(chi, 6) * Apsq[i];
                
            }

            P_vals[kk] = simpson_integrate(chi_vals, integrand);
        }
    }

    // for (const auto kRs : kRs_vals) {
    //     const auto kRs_inv = 1.0 / kRs;
    //     const auto fac2 = fac1 * power(kRs_inv, 5);
    //     const auto fac3 = beta * Rs * kRs_inv;

    //     std::vector<double> integrand(n);
    //     for (int i = 0; i < n; i++) {
    //         const auto chi = chi_vals[i];
    //         integrand[i] = fac2 * lt_dist(fac3 * chi) * power(chi, 6) * Apsq[i];
            
    //     }

    //     const auto P = simpson_integrate(chi_vals, integrand);
    //     P_vals.push_back(P);
    // }

    return PowerSpec(kRs_vals, P_vals);
}

PowerSpec zetaKin(const PowerSpec& Ekin) {
    const auto Ekin_max = Ekin.max();
    if (Ekin_max == 0.0) {
        throw std::runtime_error("Division by zero in zetaKin from Ekin.max() = 0");
    } else if (isnan(Ekin_max)) {
        throw std::runtime_error("In zetaKin: Ekin.max() = nan");
    }

    const auto zk = Ekin / Ekin_max;
    if (abs(zk.max() - 1.0) > 1e-15) {
        throw std::runtime_error("In zetaKin: Power spectrum failed normalisation test");
    }

    return zk;
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