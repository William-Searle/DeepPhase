// maths_ops.cpp
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

/*
TO DO:
- make vector class for vector arithmetic and include print_vector() function there so i can use vec.print()
*/

std::vector<double> linspace(double start, double end, std::size_t num) {
    std::vector<double> result;

    if (num == 0)
        return result;
    if (num == 1) {
        result.push_back(start);
        return result;
    }

    double step = (end - start) / (num - 1);
    result.reserve(num);

    for (std::size_t i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }

    return result;
}

std::vector<double> logspace(double start, double stop, std::size_t num) {
    std::vector<double> result;
    double log_start = std::log10(start);  // Log of the start value
    double log_stop = std::log10(stop);    // Log of the stop value

    double step = (log_stop - log_start) / (num - 1);  // Step size in log space

    for (size_t i = 0; i < num; ++i) {
        double log_val = log_start + i * step;     // Calculate log value at step i
        result.push_back(std::pow(10, log_val));   // Convert back to linear space
    }

    return result;
}

// faster than std::pow
double power(double x, int exp) {
    if (exp == 0) return 1.0;  // x^0 = 1
    double result = x;
    for (int i = 1; i < std::abs(exp); ++i) {
        result *= x;
    }
    return (exp > 0) ? result : 1.0 / result;  // Handle negative exponents
}

// for x^6 (Ekin slowed down a lot by std::pow)
double power6(double x) {
    double x2 = x * x;
    double x4 = x2 * x2;
    return x4 * x2;
}

std::string to_string_with_precision(double value, int precision) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

// simpson with non-uniform spacing - might be unstable?
double simpson_integrate(const std::vector<double>& x, const std::vector<double>& y) {
    const size_t n = x.size();
    if (n != y.size()) {
        throw std::invalid_argument("x and y must be the same size");
    }
    if (n < 2) {
        throw std::invalid_argument("Need at least two points for integration");
    }

    double integral = 0.0;
    size_t limit = (n % 2 == 0) ? n - 1 : n;  // If even, stop at n-1 to apply Simpson's

    for (size_t i = 0; i + 2 < limit; i += 2) {
        double h = x[i + 2] - x[i];
        double h1 = x[i + 1] - x[i];
        double h2 = x[i + 2] - x[i + 1];

        // If spacing is non-uniform, adjust accordingly
        if (std::abs(h1 - h2) > 1e-8) {
            // fallback to composite trapezoid for irregular spacing
            integral += 0.5 * (x[i + 1] - x[i]) * (y[i] + y[i + 1]);
            integral += 0.5 * (x[i + 2] - x[i + 1]) * (y[i + 1] + y[i + 2]);
        } else {
            integral += (h / 6.0) * (y[i] + 4 * y[i + 1] + y[i + 2]);
        }
    }

    // Trapezoid on the last interval if n is even
    if (n % 2 == 0) {
        double h = x[n - 1] - x[n - 2];
        integral += 0.5 * h * (y[n - 2] + y[n - 1]);
    }

    return integral;
}