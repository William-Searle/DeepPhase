// profile.cpp
#include <cmath>
#include <boost/numeric/odeint.hpp>
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
*/

// finite size stuff
// for (size_t i = 0; i < xi_vals.size(); ++i) {
//     assert(std::isfinite(v_vals[i]));
//     assert(std::isfinite(w_vals[i]));
//     // assert(v_vals[i] >= 0.0 && v_vals[i] <= 1.0);
//     // assert(w_vals[i] >= 0.0);  // domain-specific range?
// }

namespace plt = matplotlibcpp;
using namespace boost::numeric::odeint;

namespace Hydrodynamics { // calculate bubble profile

/*************************** Fluid profile ODE **********************************/
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

/***** FluidProfile class *****/
// Define ctor
FluidProfile::FluidProfile(const PhaseTransition::PTParams& params)
    : y0_(),
      params_(params),
      csq_(params_.csq()),
      xi_vals_(), v_vals_(), w_vals_(),
      v_prof_(), w_prof_(), la_prof_() 
    {
        // define initial state vector (xi0, v0, w0) = (vw+dlt, v+, w+)
        // is this just for bag model?
        // (starts integration just outside of wall)
        const auto dlt = 0.01;
        const auto xi0 = params_.vw() + dlt;
        y0_.push_back(xi0);

        const auto vp = params_.vpm()[0];
        const auto v0 = params_.vUF(vp);
        y0_.push_back(v0);

        const auto wp = params_.wpm()[0];
        const auto w0 = 0.1; // PLACEHOLDER
        y0_.push_back(w0);

        // define bubble profile interpolating functions
        profile();
    }

// Public functions
void FluidProfile::write(const std::string& filename) const {
    std::cout << "Writing fluid profile to disk... ";
    std::ofstream file(filename);
    file << "xi,v,w,la\n";
    // which one to use?? prof used in calculations, so maybe this?
    for (size_t i = 0; i < xi_vals_.size(); ++i) {
        file << xi_vals_[i] << "," << v_vals_[i] << "," << w_vals_[i] << "," << la_vals_[i] << "\n";
    }
    file.close();
    std::cout << "Saved to " << filename << "!\n";

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
    plt::ylim(0.0, 1.0);
    plt::grid(true);

    // w(xi)
    plt::subplot2grid(1, 3, 0, 1);
    plt::plot(xi_vals_, w_vals_);
    plt::xlabel("xi");
    plt::ylabel("w(xi)");
    plt::xlim(0.0, 1.0);
    // plt::ylim(0.0, 1.0);
    plt::grid(true);

    // la(xi)
    plt::subplot2grid(1, 3, 0, 2);
    plt::plot(xi_vals_, la_vals_, {{"linestyle", "-"}, {"color", "red"}});
    plt::plot(xi_vals_, la_vals_test_, {{"linestyle", "--"}, {"color", "blue"}});
    plt::xlabel("xi");
    plt::ylabel("la(xi)");
    plt::xlim(0.0, 1.0);
    // plt::ylim(0.0, 1.0);
    plt::grid(true);

    plt::suptitle("vw = " + to_string_with_precision(params_.vw()) + ", alpha = " + to_string_with_precision(params_.alpha()));
    plt::save(filename);

    std::cout << "Bubble profile plot saved to '" << filename << "'." << std::endl;

    return;
}

// this function might be a bit circular since it defines a FluidSystem and xi_vals, v_vals itself
// Warning: doesn't work for w profile yet!
void FluidProfile::generate_streamplot_data(int xi_pts, int y_pts, const std::string& filename) const {
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

            // double dxi = std::numeric_limits<double>::quiet_NaN();
            // double dv = std::numeric_limits<double>::quiet_NaN();

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
std::vector<state_type> FluidProfile::read(const std::string& filename) const {
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

std::vector<state_type> FluidProfile::solve_profile(int n) const {
    FluidSystem fluid(params_);
    auto y = y0_; // integrator needs non-const initial state
    
    std::cout << "(xi0,w0,w0)=(";
    for (double vals : y) {
        std::cout << vals << ",";
    }
    std::cout << ")\n";

    std::vector<state_type> states;
    std::vector<double> times;

    const auto ti = 0.0;
    const auto tf = 5.0;
    const auto dt = (tf - ti) / n;

    // integrate_const(runge_kutta4<state_type>(), fluid, y0, ti, tf, dt, push_back_state(states, times));

    runge_kutta4<state_type> stepper; // integration type
    int count = 1; // counts integration steps

    for (auto t = ti; t < tf; t += dt) {
        // std::cout << "t=" << t << ",\t" << "y=(" << y[0] << "," << y[1] << "," << y[2] << ")\n";

        states.push_back(y);
        times.push_back(t);

        const auto y_prev = y;
        stepper.do_step(fluid, y, t, dt);

        if (y[0] > 1.0 || y[0] < 0.0) {
            std::cout << "Stopping integration prematurely (n=" << count << "/" << n << ") at tau=" << t << ", xi=" << y_prev[0] << ". Further integration will yield domain error (xi < 0 or xi > 1)" << std::endl;
            break;
        }
        count++;
    }

    

    return states; // need output of times anywhere?
}

void FluidProfile::profile(bool read_prof) { // stores solve_profile vals
    if (!xi_vals_.empty()) {
        std::cerr << "Warning: xi=r/t vector not empty. Clearing values for bubble profile calculation." << std::endl;
        xi_vals_.clear();
    }

    state_type v_vals, w_vals, la_vals;

    if (read_prof) { // read fluid profile from file
        // WARNING: doesn't compare vw, alpha used in input file to params_
        const std::vector<state_type> data = read("input_profile.csv");
        xi_vals_ = data[0];
        v_vals = data[1]; // only store interpolating func values
        w_vals = data[2];
        la_vals = data[3];
    } else { // calculate fluid profile using ODE solver
        std::cout << "Warning: Profile solver not finished! Use pre-calculated bubble profile instead." << std::endl;
        
        const auto prof = solve_profile();
        for (size_t i = 0; i < prof.size(); i++) {
            xi_vals_.push_back(prof[i][0]);
            v_vals.push_back(prof[i][1]); // only store interpolating func values
            w_vals.push_back(prof[i][2]);
        }
        la_vals = calc_lambda_vals(w_vals);
    }

    // build interpolating functions
    if (v_prof_.is_initialised() || w_prof_.is_initialised() || la_prof_.is_initialised()) {
        std::cerr << "Warning: Overwriting existing interpolating functions in FluidProfile." << std::endl;
    }
    v_prof_.build(xi_vals_, v_vals);
    w_prof_.build(xi_vals_, w_vals);
    la_prof_.build(xi_vals_, la_vals);
    // v_prof_ = CubicSpline<double>(xi_vals_, v_vals);
    // w_prof_ = CubicSpline<double>(xi_vals_, w_vals);
    // la_prof_ = CubicSpline<double>(xi_vals_, la_vals);

    // store interpolated profile vals
    if (!v_vals_.empty() || !w_vals_.empty() || !la_vals_.empty()) {
        std::cerr << "Warning: Overwriting existing profile values in FluidProfile." << std::endl;
        v_vals_.clear();
        w_vals_.clear();
        la_vals_.clear();
    }

    for (const auto xi : xi_vals_) {
        v_vals_.push_back(v_prof_(xi));
        w_vals_.push_back(w_prof_(xi));
        la_vals_.push_back(la_prof_(xi));
    }

    return;
}

state_type FluidProfile::calc_lambda_vals(state_type w_vals) const { // change for when state_type = vec<double>
    if (params_.model() == "bag") {
        state_type result;
        for (const auto w : w_vals) {
            result.push_back((3.0 / 4.0) * (w - 1.0));
        }
        return result;
    } else {
        throw std::invalid_argument("In lambda(xi): only bag model implemented so far");
    }
}

// unused
interp_type FluidProfile::calc_lambda_prof() const {
    if (params_.model() == "bag") {
        return (3.0 / 4.0) * (w_prof_ - 1.0);
    } else {
        throw std::invalid_argument("In lambda(xi): only bag model implemented so far");
    }
}
/******************************/

} // namespace Hydrodynamics