// hydrodynamics.cpp

#include <cmath>
#include <vector>
#include <stdexcept>

// Speed of sound in plasma
constexpr double cs2 = 1.0 / 3.0;

// Compute mu(v) = (v - cs^2) / (1 - v^2)
inline double mu(double v) {
    return (v - cs2) / (1.0 - v * v);
}

// Compute w_out(w) = (1 + cs^2 * mu(v)) / (mu(v) + cs^2)
inline double getwow(double v) {
    double m = mu(v);
    return (1.0 + cs2 * m) / (m + cs2);
}

// Derivative dw/dv from Python dfdv
inline double dfdv(double v, double w) {
    double m = mu(v);
    double dm_dv = (1.0 - v * v + 2.0 * v * (v - cs2)) / ((1.0 - v * v) * (1.0 - v * v));
    return (w / v) * (dm_dv / (m + cs2) - (cs2 * dm_dv) / (1 + cs2 * m));
}

// Basic Runge-Kutta 4 ODE solver for dw/dv = f(v, w)
inline void solveODE(double v0, double w0, double v_end, std::vector<double>& v_values, std::vector<double>& w_values, int steps = 1000) {
    double dv = (v_end - v0) / steps;
    double v = v0;
    double w = w0;

    v_values.clear();
    w_values.clear();

    v_values.push_back(v);
    w_values.push_back(w);

    for (int i = 0; i < steps; ++i) {
        double k1 = dfdv(v, w);
        double k2 = dfdv(v + 0.5 * dv, w + 0.5 * dv * k1);
        double k3 = dfdv(v + 0.5 * dv, w + 0.5 * dv * k2);
        double k4 = dfdv(v + dv, w + dv * k3);

        w += dv * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        v += dv;

        v_values.push_back(v);
        w_values.push_back(w);
    }
}
