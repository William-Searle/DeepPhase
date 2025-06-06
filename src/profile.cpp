// profile.cpp
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>

#include "profile.hpp"
#include "PhaseTransition.hpp"
#include "hydrodynamics.hpp"
#include "physics.hpp"
#include "matplotlibcpp.h"
#include "maths_ops.hpp"

/*
TO DO:
- update generate_streamplot_data() to include points of interest (fixed pts, detonation/deflag/hybrid regions)
- fix calc of w profile in generate_streamplot_data()
- how to choose times to integrate over in solve_prof()? endpoints of integration hardcoded currently
- how to call csq? directly from PTParams or create local copy in Profile class?
- in profile, need to check for okay initial conditions otherwise v_prof etc could be divergent/solution doesn't exist
    - make timesteps dynamic if it has to stop integrating prematurely in profile_solver() (so enough integration points)
        - if it stops prematurely, find tau_max and redo the integration? this is simplest but longest solution
    - could get rid of boost integration and do the same thing generate_streamplot() does
        - i.e. start with y0 and step forwards by calculating derivatives
        - much simpler, i know it works properly (maybe takes longer, not sure?)
    - could also just specify dt rather than number of steps
    - change to integrate over specific range of xi (similar to xiao's code)
- modify plot() so you can individually plot v,w,la too (i.e. separate plot functions)
- check xi, v, w, la values should be finite and in valid range
- make xi,w,w,la vals and profiles constant after calling ctor somehow?
- move y0 definintion in ctor into a function
- might be a mistake in lambda profile calc - value for xi<vw different to xiao's code (otherwise okay)
- add 'dev mode' option for profile() where it checks lambda calculation is correct and if read profile matches ode solver profile
- profile splines don't always need to be built, only when plotting them. Get rid of build spline in ctor?
*/

// finite size stuff
// for (size_t i = 0; i < xi_vals.size(); ++i) {
//     assert(std::isfinite(v_vals[i]));
//     assert(std::isfinite(w_vals[i]));
//     // assert(v_vals[i] >= 0.0 && v_vals[i] <= 1.0);
//     // assert(w_vals[i] >= 0.0);  // domain-specific range?
// }

namespace plt = matplotlibcpp;

namespace Hydrodynamics { // calculate bubble profile

/*************************** Fluid profile ODE **********************************/
/*
EoM for a perfect fluid comes from \partial_{\mu} T^{mu nu} = 0:
    2v/xi = gamma^2 (1 - v*xi)(mu^2/cs^2 - 1) v'
    w'/w = gamma^2 * mu (1 + 1/cs^2) v'
where xi=r/t, gamma = 1/sqrt(1-v^2), v'=dv/dxi, w'=dw/dxi and mu=(xi-v)/(1-xi*v).
Parametrising these eq's gives us a set of parametric eq's (dxi/dtau,
dv/dtau and dw/dtau) to solve.
*/
double FluidSystem::dxi_dtau(double xi, double v) const {
    return xi * ((xi-v)*(xi-v) - params_.csq() * (1-xi*v)*(1-xi*v));
}

double FluidSystem::dv_dtau(double xi, double v) const {
    return 2.0 * v * params_.csq() * (1-v*v) * (1 - xi*v);
}

double FluidSystem::dw_dtau(double xi, double v, double w) const {
    return w * (1 + 1/params_.csq()) * Physics::gammaSq(v) * mu(xi, v) * dv_dtau(xi, v);
}
/*******************************************************************************/

void FluidSystem::operator()(const state_type &y, state_type &dydt, double /*tau*/) const {
    const double xi = y[0];
    const double v = y[1];
    const double w = y[2];

    dydt[0] = dxi_dtau(xi, v);
    dydt[1] = dv_dtau(xi, w);
    dydt[2] = dw_dtau(xi, v, w);

    return;
}

void push_back_state::operator()(const state_type &y, double t) const {
    // for testing
    std::cout << "t=" << t << ",\t" << "y=(";
    for (auto yy : y) {
        std::cout << yy << ",";
    }
    std::cout << ")\n";

    if (y[0] < 0.0 || y[0] >= 1.0)
        throw std::invalid_argument("Domain error in fluid profile solver: 0 <= xi < 1");
    states.push_back(y);
    times.push_back(t);
    return;
}

/****************************** FluidProfile class ******************************/
// Define ctor
FluidProfile::FluidProfile(const PhaseTransition::PTParams& params, const bool test)
    : params_(params),
      csq_(params_.csq()),
      y0_(),
      xi_vals_(), v_vals_(), w_vals_(), la_vals(),
      mode_(),
      vp_(), vpUF_(), vm_(), vmUF_(),
      wp_(), wm_()
    {
        // define hydrodynamic mode and initial state vector (xi0, v0, w0)
        get_mode_y0();

        // calculate fluid profiles v(xi), w(xi), la(xi)
        if (test) { // read fluid profile from file (testing only)
            // WARNING: doesn't compare vw, alpha used in input file to params_
            const std::vector<state_type> data = read("input_profile.csv");
            xi_vals_ = data[0];
            v_vals_ = data[1];
            w_vals_ = data[2];
            la_vals_ = data[3];
        } else { // calculate fluid profile using ODE solver
            std::cout << "Warning: Profile solver not finished! Use pre-calculated bubble profile instead." << std::endl;
            
            const auto prof = solve_profile();
            for (size_t i = 0; i < prof.size(); i++) {
                xi_vals_.push_back(prof[i][0]);
                v_vals_.push_back(prof[i][1]);
                w_vals_.push_back(prof[i][2]);
            }
            la_vals_ = calc_lambda_vals(w_vals_);
        }

        // build interpolating functions (probably not needed)
        // if (v_prof_.is_initialised() || w_prof_.is_initialised() || la_prof_.is_initialised()) {
        //     std::cerr << "Warning: Overwriting existing interpolating functions in FluidProfile." << std::endl;
        // }
        // v_prof_.build(xi_vals_, v_vals_);
        // w_prof_.build(xi_vals_, w_vals_);
        // la_prof_.build(xi_vals_, la_vals_);
    }

// Public functions
void FluidProfile::write(const std::string& filename) const {
    std::cout << "Writing fluid profile to disk... ";
    std::ofstream file(filename);
    file << "xi,v,w,la\n";

    for (size_t i = 0; i < xi_vals_.size(); ++i) {
        file << xi_vals_[i] << "," << v_vals_[i] << "," << w_vals_[i] << "," << la_vals_[i] << "\n";
    }
    file.close();

    std::cout << "Fluid profile saved to " << filename << "!\n";

    return;
}

void FluidProfile::plot(const std::string& filename) const {
    namespace plt = matplotlibcpp;

    plt::figure_size(2400, 600);

    // v(xi)
    plt::subplot2grid(1, 3, 0, 0);
    plt::plot(xi_vals_, v_vals_);
    plt::xlabel("xi");
    plt::ylabel("v(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    // w(xi)
    plt::subplot2grid(1, 3, 0, 1);
    plt::plot(xi_vals_, w_vals_);
    plt::xlabel("xi");
    plt::ylabel("w(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    // la(xi)
    plt::subplot2grid(1, 3, 0, 2);
    plt::plot(xi_vals_, la_vals_);
    plt::xlabel("xi");
    plt::ylabel("la(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    plt::suptitle("vw = " + to_string_with_precision(params_.vw()) + ", alpha = " + to_string_with_precision(params_.alpha()));
    plt::save(filename);

    std::cout << "Bubble profile plot saved to '" << filename << "'." << std::endl;

    return;
}

// Warning: doesn't work for w profile yet!
// Add lines to distinguish physical regions, shocks and separate det/deflag solutions
// This is independent of calculated profile and just needs to pass in a PTParams object. Move outside of class?
void FluidProfile::generate_streamplot_data(int xi_pts, int y_pts, const std::string& filename) const {
    std::cout << "Generating streamplot data for fluid profile... ";
    
    std::ofstream file(filename);
    file << "xi,v,w,dxidtau,dvdtau,dwdtau\n";
    file << std::fixed << std::setprecision(8); // needed for compatibility with python streamplot

    FluidSystem fluid(params_);

    // Define grid ranges (avoid's singularity at xi=0)
    const double xi_min = 0.01;
    const double xi_max = 0.99;
    const double y_min = 0.01; // bounds for v, w the same
    const double y_max = 0.99;

    // create grid for streamplot
    const auto xi_vals = linspace(xi_min, xi_max, xi_pts);
    const auto y_vals = linspace(y_min, y_max, y_pts);

    for (double xi : xi_vals) {
        for (double y : y_vals) {
            // Avoid division by zero
            if (std::abs(1 - xi * y) < 1e-6) continue;

            const auto dxi = fluid.dxi_dtau(xi, y);
            const auto dv = fluid.dv_dtau(xi, y);
            const auto dw = dv; // UPDATE THIS -> need dw at fixed v, but what v?

            file << xi << "," << y << "," << y << "," << dxi << "," << dv << "," << dw << "\n";
        }
    }

    file.close();
    std::cout << "Streamplot data saved to " << filename << "\n";
    
    return;
}

// vsh(xi) = (3 xi^2 - 1) / (2 xi) for Bag
double FluidProfile::v_shock(double xi) const {
    return (xi * xi - csq_) / ((1.0 - csq_) * xi);
}

// wsh(xi) = (9 xi - 1) / (3 (1 - xi^2)) for Bag
double FluidProfile::w_shock(double xi) const {
    return (xi - csq_ * csq_) / (csq_ * (1.0 - xi * xi));
}

// Private functions

// determine hydrodynamic mode
void FluidProfile::get_mode_y0() const {
    const auto vw = params_.vw();
    const auto vwsq = vw * vw;

    const auto cmsq = params_.cmsq();
    const auto cpsq = params_.cpsq();

    // if (vwsq < cpsq) { // deflagration
    //     mode_ = 0;
    //     y0_[0] = xi_shock(); // xi0 - start integration at shock

    //     vmUF_ = 0.0;
    //     vm_ = vw;
    //     vp = ; // need to be careful calculating this since alpha /= alpha+
    //     vpUF_ = mu(vw, vp_);
    //     y0_[1] = v_shock(vpUF_); // v0 = v(xi_sh) = v1 (fluid velocity just behind shock)

    // }
    const auto dlt = 0.01; // wall and shock are discontinuities for v(xi) so start integration just before them
    if (vwsq >= cpsq) { // detonation
        mode_ = 2;
        y0_[0] = vw - dlt; // xi0 = xi_w

        vpUF_ = 0.0; // v+ (centre of bubble frame / universe frame)
        vp_ = vw; // v+ (bubble wall frame)
        vm_ = calc_vm(vp_);
        vmUF_ = mu(vw, vm_);
        y0_[1] = vmUF_; // v0 = v(xi_w) = vm(UF)

        // wp_ = params_.wN();
        wp_ = 1.0;
        wm_ = calc_wm(wp_, vp_, vm_);
        y0_[2] = wm_; // w0 = w(xi_w) = wm
    } else {
        throw std::invalid_argument("Deflagration/Hybrid solutions not yet implemented");
    }


}

// sgn=1 for detonation, not sure what condition specifically fixes this though
// sgn=-1 for deflag/hybrid?
// this is for bag model only, generalise for any EoS?
double FluidProfile::calc_vm(double vp) const { // vm from vp
    const auto sgn = 1.0;
    const auto fac = vp * (1.0 + alpha_) / 2.0 + (1./3. - alpha_) / (2.0 * vp);
    return fac + sgn * std::sqrt(fac * fac - 1./3.);
}

// bag model only
double FluidProfile::calc_vp(double vm) const { // vp from vm
    const auto sgn = 1.0;
    const auto fac = vm / 2.0 + 1.0 / (6.0 * vm);
    return (fac + sgn * std::sqrt(fac * fac + alpha_ * alpha_ + (2./3.) * alpha_ - 1./3.)) / (1.0 + alpha_);
}

double FluidProfile::calc_wm(double wp, double vp, double vm) const {
    // from matching condition -> generic for all EoS and hydrodynamic models using perfect fluid EM tensor
    return wp * vp * (1.0 - vm * vm) / (vm * (1.0 - vp * vp));
}

std::vector<state_type> FluidProfile::read(const std::string& filename) const {
    std::cout << "Warning: Read fluid profile does not check PT parameters of input file. Manual entry of PT parameters required!\n";

    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open file " + filename);
    }

    std::string line;
    std::getline(file, line); // Skip header

    state_type xi_vals, v_vals, w_vals, la_vals;

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::array<double, 4> values;
        std::string token;

        for (auto& val : values) {
            if (!std::getline(ss, token, ',')) {
                throw std::runtime_error("Malformed line in " + filename + ": " + line);
            }
            val = std::stod(token);
        }

        xi_vals.push_back(values[0]);
        v_vals.push_back(values[1]);
        w_vals.push_back(values[2]);
        la_vals.push_back(values[3]);
    }

    return {xi_vals, v_vals, w_vals, la_vals};
}

// This is for Bag EoS
state_type FluidProfile::calc_lambda_vals(state_type w_vals) const {
    state_type result;
    for (const auto w : w_vals) {
        result.push_back((3.0 / 4.0) * (w - 1.0));
    }
    return result;
}

std::vector<state_type> FluidProfile::solve_profile(int n) const {    
    std::cout << "Solving fluid profile for hydrodynamic mode=";
    if (mode_ == 0) {
        std::cout << "deflagration...\n";
    } else if (mode_ == 1) {
        std::cout << "hybrid...\n";
    } else {
        std::cout << "detonation...\n";
    }

    FluidSystem fluid(params_);
    std::vector<state_type> states;
    std::vector<double> times;

    const auto ti = 0.0;
    const auto tf = 5.0;
    const auto dt = (tf - ti) / n;

    return states; // need output of times anywhere?
}
/*******************************************************************************/

} // namespace Hydrodynamics