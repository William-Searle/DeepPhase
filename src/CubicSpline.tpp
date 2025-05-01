#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>

#include "maths_ops.hpp"

// dflt ctor needed for FluidProfile class
template <typename T>
CubicSpline<T>::CubicSpline() 
    : initialized_(false) {}

template <typename T>
CubicSpline<T>::CubicSpline(const std::vector<T>& x, const std::vector<T>& y) {
    if (x.size() != y.size() || x.size() < 2)
        throw std::invalid_argument("Vectors must be the same size and have at least two points.");

    size_t n = x.size();
    this->x = x;
    this->y = y;
    h.resize(n - 1);
    a = y;
    c.resize(n);

    // Step 1: Compute h[i]
    for (size_t i = 0; i < n - 1; ++i)
        h[i] = x[i + 1] - x[i];

    // Step 2: Build the tridiagonal system
    std::vector<T> alpha(n - 1);
    for (size_t i = 1; i < n - 1; ++i)
        alpha[i] = (3 / h[i]) * (a[i + 1] - a[i]) - (3 / h[i - 1]) * (a[i] - a[i - 1]);

    std::vector<T> l(n), mu(n), z(n);
    l[0] = 1;
    mu[0] = z[0] = 0;

    for (size_t i = 1; i < n - 1; ++i) {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1;
    z[n - 1] = c[n - 1] = 0;

    b.resize(n - 1);
    d.resize(n - 1);

    for (int j = static_cast<int>(n) - 2; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }
}

template <typename T>
T CubicSpline<T>::operator()(T xi) const {
    if (!initialized_)
        throw std::runtime_error("CubicSpline was not initialized!");

    if (xi <= x.front()) return y.front();
    if (xi >= x.back()) return y.back();

    auto it = std::upper_bound(x.begin(), x.end(), xi);
    size_t i = std::distance(x.begin(), it) - 1;

    T dx = xi - x[i];
    return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
}