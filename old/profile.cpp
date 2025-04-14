// profile.cpp

#include "hydrodynamics.hpp"
#include <vector>
#include <utility>
#include <cmath>

struct Profile {
    std::vector<double> v_vals;
    std::vector<double> w_vals;
};

// Solve the ODE and return profile (w(v)) as a struct
inline Profile generate_profile(double v0, double w0, double v_end) {
    std::vector<double> v_vals;
    std::vector<double> w_vals;
    solveODE(v0, w0, v_end, v_vals, w_vals);
    return {v_vals, w_vals};
}

// Logarithmic version of the profile, for spacing similar to np.geomspace
inline Profile generate_log_profile(double v0, double w0, double v_end, int n_points = 100) {
    std::vector<double> v_vals;
    std::vector<double> w_vals;

    double log_v0 = std::log(v0);
    double log_v1 = std::log(v_end);

    for (int i = 0; i < n_points; ++i) {
        double lv = log_v0 + (log_v1 - log_v0) * i / (n_points - 1);
        v_vals.push_back(std::exp(lv));
    }

    // Assume we start from initial (v0, w0) and integrate ODE step by step
    w_vals.push_back(w0);
    for (size_t i = 1; i < v_vals.size(); ++i) {
        double dv = v_vals[i] - v_vals[i - 1];
        double v_mid = 0.5 * (v_vals[i] + v_vals[i - 1]);
        double w_prev = w_vals.back();
        double dw = dv * dfdv(v_mid, w_prev);  // simple Euler step
        w_vals.push_back(w_prev + dw);
    }

    return {v_vals, w_vals};
}
