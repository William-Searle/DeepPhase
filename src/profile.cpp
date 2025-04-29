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

// void plot_velocity_profile(const std::string& filename) {
//     std::ifstream file(filename);
//     std::string line;
//     std::vector<double> xi, v;

//     // Skip header
//     std::getline(file, line);

//     // Read values
//     while (std::getline(file, line)) {
//         std::stringstream ss(line);
//         std::string tau_str, xi_str, v_str, w_str;
//         std::getline(ss, tau_str, ',');
//         std::getline(ss, xi_str, ',');
//         std::getline(ss, v_str, ',');
//         std::getline(ss, w_str, ',');

//         xi.push_back(std::stod(xi_str));
//         v.push_back(std::stod(v_str));
//     }

//     // Plot as a vector field (xi vs v)
//     std::vector<double> u(xi.size(), 0.0); // zero horizontal component
//     std::vector<double> x = xi;
//     std::vector<double> y(xi.size(), 0.0); // all arrows start at y=0
//     std::vector<double> v_dir = v;         // vertical component

//     plt::figure_size(800, 400);
//     plt::quiver(x, y, u, v_dir, 10); // scale=10 can be adjusted
//     plt::xlabel("xi");
//     plt::ylabel("v");
//     plt::title("Velocity profile v(xi)");
//     plt::save("v_vs_xi.png");
//     std::cout << "Saved plot to v_vs_xi.png\n";
// }


} // namespace FluidProfile