#pragma once

#include <algorithm>
#include <vector>
#include <stdexcept>

template <typename T>
CubicSpline<T>::CubicSpline() : initialized_(false) {}

template <typename T>
CubicSpline<T>::CubicSpline(const std::vector<T>& x, const std::vector<T>& y) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::invalid_argument("x and y must be the same size and contain at least two points.");
    }

    const size_t n = x.size();

    // Check for strictly increasing x
    for (size_t i = 1; i < n; ++i) {
        if (x[i] <= x[i - 1]) {
            throw std::invalid_argument("x-values must be strictly increasing.");
        } else if (isnan(x[i])) {
            throw std::invalid_argument("NaN x-value detected.");
        }
    }

    x_ = x;
    y_ = y;
    h_.resize(n - 1);
    a_ = y;
    c_.resize(n);

    // Step 1: Compute h[i]
    for (size_t i = 0; i < n - 1; ++i) {
        h_[i] = x_[i + 1] - x_[i];
    }

    // Step 2: Compute alpha
    std::vector<T> alpha(n - 1, 0);
    for (size_t i = 1; i < n - 1; ++i) {
        alpha[i] = (3.0 / h_[i]) * (a_[i + 1] - a_[i]) - (3.0 / h_[i - 1]) * (a_[i] - a_[i - 1]);
    }

    // Step 3: Solve tridiagonal system
    std::vector<T> l(n), mu(n), z(n);
    l[0] = 1;
    mu[0] = z[0] = 0;

    for (size_t i = 1; i < n - 1; ++i) {
        l[i] = 2 * (x_[i + 1] - x_[i - 1]) - h_[i - 1] * mu[i - 1];
        mu[i] = h_[i] / l[i];
        z[i] = (alpha[i] - h_[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1;
    z[n - 1] = c_[n - 1] = 0;

    b_.resize(n - 1);
    d_.resize(n - 1);

    for (int j = static_cast<int>(n) - 2; j >= 0; --j) {
        c_[j] = z[j] - mu[j] * c_[j + 1];
        b_[j] = (a_[j + 1] - a_[j]) / h_[j] - h_[j] * (c_[j + 1] + 2 * c_[j]) / 3;
        d_[j] = (c_[j + 1] - c_[j]) / (3 * h_[j]);
    }

    initialized_ = true;
}

template <typename T>
T CubicSpline<T>::operator()(T xi) const {
    if (!initialized_) {
        throw std::runtime_error("CubicSpline has not been initialized.");
    }

    // Extrapolate flat for out-of-bounds input
    if (xi <= x_.front()) return y_.front();
    if (xi >= x_.back()) return y_.back();

    // Binary search to find correct interval
    auto it = std::upper_bound(x_.begin(), x_.end(), xi);
    size_t i = std::distance(x_.begin(), it) - 1;

    T dx = xi - x_[i];
    return a_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
}