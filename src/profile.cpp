// profile.cpp
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>

#include "profile.hpp"
#include "physics.hpp"
#include "matplotlibcpp.h"
#include "maths_ops.hpp"

namespace plt = matplotlibcpp;
using namespace boost::numeric::odeint;

namespace FluidProfile { // calculate bubble profile (NOT FINISHED)

/*************************** Fluid profile ODE **********************************/
double FluidSystem::dxi_dtau(double xi, double v, double csq) const {
    return xi * ((xi-v)*(xi-v) - csq * (1-xi*v)*(1-xi*v));
}

double FluidSystem::dv_dtau(double xi, double v, double csq) const {
    return 2.0 * v * csq * (1-v*v) * (1 - xi*v);
}

double FluidSystem::dw_dtau(double xi, double v, double w, double csq) const {
    const auto mu = (xi - v) / (1 - xi * v);
    return w * (1 + 1/csq) * Physics::gammaSq(v) * mu * dv_dtau(xi, v, csq);
}
/*******************************************************************************/

void FluidSystem::operator()(const state_type &y, state_type &dydt, double /*tau*/) const {
    const double xi = y[0];
    const double v = y[1];
    const double w = y[2];

    dydt[0] = dxi_dtau(xi, v, csq_);
    dydt[1] = dv_dtau(xi, w, csq_);
    dydt[2] = dw_dtau(xi, v, w, csq_);
}

void push_back_state::operator()(const state_type &y, double t) const {
    states.push_back(y);
    times.push_back(t);
}

// y0 = {xi0, v0, w0} (vector of initial values)
void solve_prof(state_type y0, double csq, int n) {
    FluidSystem fluid(csq);

    std::vector<state_type> states;
    std::vector<double> times;

    const auto ti = 0.0;
    const auto tf = 5.0;
    const auto dt = (tf - ti) / n;

    integrate_const(runge_kutta4<state_type>(), fluid, y0, ti, tf, dt, push_back_state(states, times));

    // Save to CSV
    std::ofstream file("fluid_solution.csv");
    file << "tau,xi,v,w\n";
    for (size_t i = 0; i < times.size(); ++i) {
        file << times[i] << "," << states[i][0] << "," << states[i][1] << "," << states[i][2] << "\n";
    }
    file.close();
    std::cout << "Saved to fluid_solution.csv\n";

    return;
}

void generate_streamplot_data(double csq, int xi_pts, int v_pts, const std::string& filename) {
    std::ofstream file(filename);
    file << "xi,v,dxidtau,dvdtau\n";
    file << std::fixed << std::setprecision(8); // needed for compatibility with python streamplot

    FluidProfile::FluidSystem fluid(csq);

    // Define the ranges for xi and v (avoid's singularity at xi=0)
    const double xi_min = 0.01;
    const double xi_max = 0.99;
    const double v_min = 0.01;
    const double v_max = 0.99;

    // create grid for streamplot
    const auto xi_vals = linspace(xi_min, xi_max, xi_pts);
    const auto v_vals = linspace(v_min, v_max, v_pts);

    for (double xi : xi_vals) {
        for (double v : v_vals) {

            // double dxi = std::numeric_limits<double>::quiet_NaN();
            // double dv = std::numeric_limits<double>::quiet_NaN();

            // Avoid division by zero
            if (std::abs(1 - xi * v) < 1e-6) continue;

            const auto dxi = fluid.dxi_dtau(xi, v, csq);
            const auto dv = fluid.dv_dtau(xi, v, csq);

            file << xi << "," << v << "," << dxi << "," << dv << "\n";
        }
    }

    file.close();
    std::cout << "Streamplot data saved to " << filename << "\n";
}

} // namespace FluidProfile