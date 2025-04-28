// profile.cpp
#include <cmath>

#include "rk4_solver.hpp"
#include "physics.hpp"

namespace FluidProfile { // calculate bubble profile (NOT FINISHED)

double dxi_dtau(double xi, double v, double csq) {
    return xi * (std::pow(xi-v,2) - csq * std::pow(1-xi*v,2));
}

double dv_dtau(double xi, double v, double csq) {
    return 2.0 * v * csq * (1-std::pow(v,2)) * (1 - xi*v);
}

double dw_dtau(double xi, double v, double w, double csq) {
    const auto mu = (xi - v) / (1 - xi * v);
    const auto dvdtau = dv_dtau(xi, v, csq);
    return w * (1 + 1/csq) * Physics::gammaSq(v) * mu * dvdtau;
}

RK4::State profile(RK4::State& init, double t0, double tf, int n) {
    RK4::DerivFunc system = [](const RK4::State& s, double /*tau*/) {
        double xi = s.xi;
        double v = s.v;
        double w = s.w;

        double csq = 1.0/3.0; // change to dflt value in lambda func

        double dxidtau = dxi_dtau(xi, v, csq);
        double dvdtau = dv_dtau(xi, v, csq);
        double dwdtau = dw_dtau(xi, v, w, csq);

        return RK4::State{dxidtau, dvdtau, dwdtau};
    };

    double dt = (tf - t0) / n;

    return RK4::integrate(system, init, t0, tf, dt);

}

} // namespace FluidProfile