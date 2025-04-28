// rk4_solver.cpp
#include "rk4_solver.hpp"

/*
TO DO:
- change xi, v, w in structure to x1, x2, x3
*/

namespace RK4 {

State State::operator+(const State& other) const {
    return { xi + other.xi, v + other.v, w + other.w };
}

State State::operator*(double scalar) const {
    return { scalar * xi, scalar * v, scalar * w };
}

State integrate(const DerivFunc& f, const State& v0, double t0, double tf, double dt) {
    State v = v0;
    double t = t0;

    while (t < tf) {
        State k1 = f(v, t);
        State k2 = f(v + k1 * (dt / 2.0), t + dt / 2.0);
        State k3 = f(v + k2 * (dt / 2.0), t + dt / 2.0);
        State k4 = f(v + k3 * dt, t + dt);

        v = v + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
        t += dt;
    }

    return v;
}

} // namespace RK4
