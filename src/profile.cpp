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
- define lambda in profile() (currently has placeholder so integration routine can be tested)
- update generate_streamplot_data() to include points of interest (fixed pts, detonation/deflag/hybrid regions)
- fix calc of w profile in generate_streamplot_data()
- how to choose times to integrate over in solve_prof()? endpoints of integration hardcoded currently
- how to call csq? directly from PTParams or create local copy in Profile class?
- in profile, need to check for okay initial conditions otherwise v_prof etc could be divergent/solution doesn't exist
- make timesteps dynamic if it has to stop integrating prematurely in profile_solver() (so enough integration points)
    - if it stops prematurely, find tau_max and redo the integration? this is simplest but longest solution
*/

namespace plt = matplotlibcpp;
using namespace boost::numeric::odeint;

namespace Hydrodynamics { // calculate bubble profile

/*************************** Fluid profile ODE **********************************/
double FluidSystem::dxi_dtau(double xi, double v) const {
    return xi * ((xi-v)*(xi-v) - csq_ * (1-xi*v)*(1-xi*v));
}

double FluidSystem::dv_dtau(double xi, double v) const {
    return 2.0 * v * csq_ * (1-v*v) * (1 - xi*v);
}

double FluidSystem::dw_dtau(double xi, double v, double w) const {
    return w * (1 + 1/csq_) * Physics::gammaSq(v) * mu(xi, v) * dv_dtau(xi, v);
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
FluidProfile::FluidProfile(PhaseTransition::PTParams& params)
    : y0_(),
      csq_( params.cmsq() ), // bad to do this? just call params.cmsq() when needed instead?
      xi_vals_(), v_vals_(), w_vals_(), la_vals_(),
      v_prof_(), w_prof_(), la_prof_() 
    {
        // define initial state vector (xi0, v0, w0) = (vw+dlt, v+, w+)
        // (starts integration just outside of wall)
        const auto dlt = 0.01;
        const auto xi0 = params.vw() + dlt;
        y0_.push_back(xi0);

        const auto vp = params.vpm()[0];
        const auto v0 = params.vUF(vp);
        y0_.push_back(v0);

        const auto wp = params.wpm()[0];
        const auto w0 = 0.1; // PLACEHOLDER
        y0_.push_back(w0);

        // define bubble profile interpolating functions
        // individual vals stored using profile()
        const auto prof_interp = profile();
        v_prof_ = prof_interp[0];
        w_prof_ = prof_interp[1];
        la_prof_ = prof_interp[2];
    }

// Public functions
void FluidProfile::write() const {
    std::cout << "Writing fluid profile to disk...\n";
    std::ofstream file("fluid_solution.csv");
    file << "xi,v,w\n";
    for (size_t i = 0; i < xi_vals_.size(); ++i) {
        file << xi_vals_[i] << "," << v_vals_[i] << "," << w_vals_[i] << "\n";
    }
    file.close();
    std::cout << "Saved to fluid_solution.csv\n";

    return;
}

void FluidProfile::plot(const std::string& filename) const {
    namespace plt = matplotlibcpp;

    plt::figure_size(1600, 600);

    // First plot (top half)
    plt::subplot2grid(1, 2, 0, 0);
    plt::plot(xi_vals_, v_vals_);
    plt::xlabel("xi");
    plt::ylabel("v(xi)");
    plt::xlim(0.0, 1.0);
    plt::ylim(0.0, 1.0);
    plt::grid(true);

    // Second plot (bottom half)
    plt::subplot2grid(1, 2, 0, 1);
    plt::plot(xi_vals_, w_vals_);
    plt::xlabel("xi");
    plt::ylabel("w(xi)");
    plt::xlim(0.0, 1.0);
    plt::ylim(0.0, 1.0);
    plt::grid(true);

    // suptitle("(vw, alpha) = " + );
    plt::save(filename);

    return;
}

// this function might be a bit circular since it defines a FluidSystem and xi_vals, v_vals itself
// Warning: doesn't work for w profile yet!
void FluidProfile::generate_streamplot_data(int xi_pts, int y_pts, const std::string& filename) const {
    std::ofstream file(filename);
    file << "xi,v,w,dxidtau,dvdtau,dwdtau\n";
    file << std::fixed << std::setprecision(8); // needed for compatibility with python streamplot

    FluidSystem fluid(csq_);

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

// Private functions
std::vector<state_type> FluidProfile::solve_profile(int n) const {
    FluidSystem fluid(csq_);
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

std::vector<CubicSpline<double>> FluidProfile::profile() {
    const auto prof = solve_profile();
    std::vector<double> xi_vals, v_vals, w_vals, la_vals;

    std::ofstream file("test_solver.csv");
    file << "xi,v,w\n";

    // fill vectors for xi, v, w
    for (size_t i = 0; i < prof.size(); i++) {
        xi_vals.push_back(prof[i][0]);
        v_vals.push_back(prof[i][1]);
        w_vals.push_back(prof[i][2]);

        file << xi_vals[i] << "," << v_vals[i] << "," << w_vals[i] << "\n";
    }
    file.close();

    // define lambda here, placeholder for now!!
    // make separate function to calc lambda?
    // interp for both la and w or define la_prof directly from w_prof?
    la_vals = v_vals;

    // store xi, v, w values in class
    // not sure if this is the best implementation
    xi_vals_ = xi_vals;
    v_vals_ = v_vals;
    w_vals_ = w_vals;
    la_vals_ = la_vals;

    // define interpolating functions
    CubicSpline<double> v_prof(xi_vals, v_vals);
    CubicSpline<double> w_prof(xi_vals, w_vals);
    CubicSpline<double> la_prof(xi_vals, la_vals);

    // vector of interpolating functions
    std::vector<CubicSpline<double>> prof_interp;
    prof_interp.push_back(v_prof);
    prof_interp.push_back(w_prof);
    prof_interp.push_back(la_prof);

    return prof_interp; // {v(xi), w(xi), la(xi)}
}
/******************************/

} // namespace Hydrodynamics