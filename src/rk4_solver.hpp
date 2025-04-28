// rk4_solver.hpp
#pragma once

#include <functional>

namespace RK4 {

struct State {
    double xi;
    double v;
    double w;

    State operator+(const State& other) const;
    State operator*(double scalar) const;
};

using DerivFunc = std::function<State(const State&, double)>;

State integrate(const DerivFunc& f, const State& y0, double t0, double tf, double dt);

} // namespace RK4