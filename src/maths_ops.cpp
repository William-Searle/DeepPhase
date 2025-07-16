// maths_ops.cpp
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <functional>
#include <cassert>
// #include <matplotlibcpp.h>

#include "constants.hpp"
#include "maths_ops.hpp"

/*
TO DO:
- make vector class for vector arithmetic and include print_vector() function there so i can use vec.print()
*/

// void plot(std::vector<double> x_vals, std::vector<double> y_vals, const std::string& filename) {
//     namespace plt = matplotlibcpp;

//     plt::figure_size(800, 600);
//     plt::plot(x_vals, y_vals);
//     plt::grid(true);
//     plt::save("../" + filename);

//     return;
// }

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

double power3(double x) {
    return x * x * x;
}

std::string to_string_with_precision(double value, int precision) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

// simpson with non-uniform spacing - might be unstable?
// add error estimation somehow... (or maybe check against boost.math integration)
// could write check_stability() function that performs the integration a number of times with different number of pts to check when it converges
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

double simpson_nonuniform(const std::vector<double>& x, const std::vector<double>& y) {
    const size_t n = x.size();
    double integral = 0.0;

    size_t i = 0;
    while (i + 2 < n) {
        const double x0 = x[i],     x1 = x[i+1],     x2 = x[i+2];
        const double y0 = y[i],     y1 = y[i+1],     y2 = y[i+2];
        const double h0 = x1 - x0,  h1 = x2 - x1;

        // Use Simpson-like Lagrange formula
        const double denom = h0 * h1 * (h0 + h1);
        const double A = -h1 * h1 / denom;
        const double B = (h1 * h1 - h0 * h0) / denom;
        const double C = h0 * h0 / denom;
        integral += (A * y0 + B * y1 + C * y2) * (x2 - x0) / 2.0;

        i += 2; // advance in steps of 2
    }

    // Fallback for remaining points using trapezoid
    if (i + 1 < n) {
        integral += 0.5 * (x[i + 1] - x[i]) * (y[i] + y[i + 1]);
    }

    return integral;
}

// can't just use a wrapper for 2d vector simpson since it slows down dlt a lot
double simpson_2d_integrate_flat(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& f_flat) {
    const size_t nx = x.size();
    const size_t ny = y.size();

    if (f_flat.size() != ny * nx) {
        throw std::invalid_argument("Size of f_flat must be y.size() * x.size()");
    }

    auto f = [&](size_t j, size_t i) -> double {
        return f_flat[j * nx + i];
    };

    double total = 0.0;

    size_t nx_lim = (nx % 2 == 0) ? nx - 1 : nx;
    size_t ny_lim = (ny % 2 == 0) ? ny - 1 : ny;

    for (size_t j = 0; j + 2 < ny_lim; j += 2) {
        double hy1 = y[j + 1] - y[j];
        double hy2 = y[j + 2] - y[j + 1];
        double hy = y[j + 2] - y[j];

        for (size_t i = 0; i + 2 < nx_lim; i += 2) {
            double hx1 = x[i + 1] - x[i];
            double hx2 = x[i + 2] - x[i + 1];
            double hx = x[i + 2] - x[i];

            if (std::abs(hx1 - hx2) > 1e-8 || std::abs(hy1 - hy2) > 1e-8) {
                double area =
                    0.25 * (x[i + 1] - x[i]) * (y[j + 1] - y[j]) * (f(j, i) + f(j + 1, i) + f(j, i + 1) + f(j + 1, i + 1)) +
                    0.25 * (x[i + 2] - x[i + 1]) * (y[j + 1] - y[j]) * (f(j, i + 1) + f(j + 1, i + 1) + f(j, i + 2) + f(j + 1, i + 2)) +
                    0.25 * (x[i + 1] - x[i]) * (y[j + 2] - y[j + 1]) * (f(j + 1, i) + f(j + 2, i) + f(j + 1, i + 1) + f(j + 2, i + 1)) +
                    0.25 * (x[i + 2] - x[i + 1]) * (y[j + 2] - y[j + 1]) * (f(j + 1, i + 1) + f(j + 2, i + 1) + f(j + 1, i + 2) + f(j + 2, i + 2));
                total += area;
            } else {
                total += (hx * hy / 36.0) * (
                    f(j, i) + 4 * f(j, i + 1) + f(j, i + 2) +
                    4 * (f(j + 1, i) + 4 * f(j + 1, i + 1) + f(j + 1, i + 2)) +
                    f(j + 2, i) + 4 * f(j + 2, i + 1) + f(j + 2, i + 2)
                );
            }
        }
    }

    // Handle remaining strip in x if nx is even
    if (nx % 2 == 0) {
        for (size_t j = 0; j + 1 < ny; ++j) {
            double hy = y[j + 1] - y[j];
            double hx = x[nx - 1] - x[nx - 2];
            total += 0.25 * hx * hy * (
                f(j, nx - 2) + f(j + 1, nx - 2) +
                f(j, nx - 1) + f(j + 1, nx - 1)
            );
        }
    }

    // Handle remaining strip in y if ny is even
    if (ny % 2 == 0) {
        for (size_t i = 0; i + 1 < nx; ++i) {
            double hx = x[i + 1] - x[i];
            double hy = y[ny - 1] - y[ny - 2];
            total += 0.25 * hx * hy * (
                f(ny - 2, i) + f(ny - 2, i + 1) +
                f(ny - 1, i) + f(ny - 1, i + 1)
            );
        }
    }

    return total;
}

// unused
double simpson_2d_integrate(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f) {
    const size_t nx = x.size();
    const size_t ny = y.size();

    if (f.size() != ny || f[0].size() != nx) {
        throw std::invalid_argument("Size of f must be [y.size()][x.size()]");
    }

    double total = 0.0;

    size_t nx_lim = (nx % 2 == 0) ? nx - 1 : nx;
    size_t ny_lim = (ny % 2 == 0) ? ny - 1 : ny;

    for (size_t j = 0; j + 2 < ny_lim; j += 2) {
        double hy1 = y[j + 1] - y[j];
        double hy2 = y[j + 2] - y[j + 1];
        double hy = y[j + 2] - y[j];

        for (size_t i = 0; i + 2 < nx_lim; i += 2) {
            double hx1 = x[i + 1] - x[i];
            double hx2 = x[i + 2] - x[i + 1];
            double hx = x[i + 2] - x[i];

            if (std::abs(hx1 - hx2) > 1e-8 || std::abs(hy1 - hy2) > 1e-8) {
                // Use composite trapezoidal if spacing is non-uniform
                double area =
                0.25 * (x[i + 1] - x[i]) * (y[j + 1] - y[j]) * (f[j][i] + f[j + 1][i] + f[j][i + 1] + f[j + 1][i + 1]) +
                0.25 * (x[i + 2] - x[i + 1]) * (y[j + 1] - y[j]) * (f[j][i + 1] + f[j + 1][i + 1] + f[j][i + 2] + f[j + 1][i + 2]) +
                0.25 * (x[i + 1] - x[i]) * (y[j + 2] - y[j + 1]) * (f[j + 1][i] + f[j + 2][i] + f[j + 1][i + 1] + f[j + 2][i + 1]) +
                0.25 * (x[i + 2] - x[i + 1]) * (y[j + 2] - y[j + 1]) * (f[j + 1][i + 1] + f[j + 2][i + 1] + f[j + 1][i + 2] + f[j + 2][i + 2]);
                total += area;
            } else {
                // Use 2D Simpson's rule
                total += (hx * hy / 36.0) * (
                f[j][i] + 4 * f[j][i + 1] + f[j][i + 2] +
                4 * (f[j + 1][i] + 4 * f[j + 1][i + 1] + f[j + 1][i + 2]) +
                f[j + 2][i] + 4 * f[j + 2][i + 1] + f[j + 2][i + 2]
                );
            }
        }
    }

    // Handle remaining strip in x if nx is even
    if (nx % 2 == 0) {
        for (size_t j = 0; j + 1 < ny; ++j) {
            double hy = y[j + 1] - y[j];
            double hx = x[nx - 1] - x[nx - 2];
            total += 0.25 * hx * hy * (
            f[j][nx - 2] + f[j + 1][nx - 2] +
            f[j][nx - 1] + f[j + 1][nx - 1]);
        }
    }

    // Handle remaining strip in y if ny is even
    if (ny % 2 == 0) {
        for (size_t i = 0; i + 1 < nx; ++i) {
            double hx = x[i + 1] - x[i];
            double hy = y[ny - 1] - y[ny - 2];
            total += 0.25 * hx * hy * (
            f[ny - 2][i] + f[ny - 2][i + 1] +
            f[ny - 1][i] + f[ny - 1][i + 1]);
        }
    }

    return total;
}

// Helper to compute weights along one dimension for non-uniform points
// Inputs:
//   coords: vector of coordinates (x or y)
// Outputs:
//   weights: vector of vectors of 3 weights per interval
//   intervals: vector of interval sizes (x2 - x0)
static void precompute_1d_weights(
    const std::vector<double>& coords,
    std::vector<std::vector<double>>& weights,
    std::vector<double>& intervals)
{
    const size_t n = coords.size();
    if (n < 3) {
        throw std::invalid_argument("Grid must have at least 3 points for Simpson's rule");
    }

    weights.resize(n - 2);
    intervals.resize(n - 2);

    for (size_t i = 0; i + 2 < n; ++i) {
        const double x0 = coords[i];
        const double x1 = coords[i+1];
        const double x2 = coords[i+2];

        const double hx0 = x1 - x0;
        const double hx1 = x2 - x1;
        const double denom = hx0 * hx1 * (hx0 + hx1);

        // Compute weights for points i, i+1, i+2
        double A0 = -hx1 * hx1 / denom;
        double A1 = (hx1 * hx1 - hx0 * hx0) / denom;
        double A2 = hx0 * hx0 / denom;

        weights[i] = {A0, A1, A2};
        intervals[i] = x2 - x0;
    }
}

// define this without function call to 1d version for speed
SimpsonWeights2D precompute_simpson_weights_2d(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    SimpsonWeights2D w;

    precompute_1d_weights(x, w.Ax_weights, w.dx);
    precompute_1d_weights(y, w.Ay_weights, w.dy);

    return w;
}

// assumes log-spaced grid
double simpson_2d_nonuniform_flat(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& f_flat) {
    const size_t nx = x.size();
    const size_t ny = y.size();

    if (f_flat.size() != nx * ny) {
        throw std::invalid_argument("Size of f_flat must be x.size() * y.size()");
    }

    double total = 0.0;

    size_t i = 0;
    while (i + 2 < nx) {
        size_t j = 0;
        while (j + 2 < ny) {
            const double x0 = x[i], x1 = x[i+1], x2 = x[i+2];
            const double y0 = y[j], y1 = y[j+1], y2 = y[j+2];

            const double hx0 = x1 - x0, hx1 = x2 - x1;
            const double hy0 = y1 - y0, hy1 = y2 - y1;

            const double denom_x = hx0 * hx1 * (hx0 + hx1);
            const double Ax0 = -hx1 * hx1 / denom_x;
            const double Ax1 = (hx1 * hx1 - hx0 * hx0) / denom_x;
            const double Ax2 = hx0 * hx0 / denom_x;
            const double dx = x2 - x0;

            const double denom_y = hy0 * hy1 * (hy0 + hy1);
            const double Ay0 = -hy1 * hy1 / denom_y;
            const double Ay1 = (hy1 * hy1 - hy0 * hy0) / denom_y;
            const double Ay2 = hy0 * hy0 / denom_y;
            const double dy = y2 - y0;

            const size_t row0 = j * nx;
            const size_t row1 = (j+1) * nx;
            const size_t row2 = (j+2) * nx;

            double local = 0.0;
            local += Ax0 * Ay0 * f_flat[row0 + i];
            local += Ax1 * Ay0 * f_flat[row0 + i+1];
            local += Ax2 * Ay0 * f_flat[row0 + i+2];

            local += Ax0 * Ay1 * f_flat[row1 + i];
            local += Ax1 * Ay1 * f_flat[row1 + i+1];
            local += Ax2 * Ay1 * f_flat[row1 + i+2];

            local += Ax0 * Ay2 * f_flat[row2 + i];
            local += Ax1 * Ay2 * f_flat[row2 + i+1];
            local += Ax2 * Ay2 * f_flat[row2 + i+2];

            total += dx * dy * local / 4.0;
            j += 2;
        }

        // Trapezoidal rule for last row if ny is even
        if (ny % 2 == 0) {
            const double hy = y[ny - 1] - y[ny - 2];
            const size_t row0 = (ny - 2) * nx;
            const size_t row1 = (ny - 1) * nx;
            for (size_t ii = 0; ii + 1 < 3; ++ii) {
                const double hx = x[i + ii + 1] - x[i + ii];
                total += 0.25 * hx * hy * (
                    f_flat[row0 + i + ii] + f_flat[row1 + i + ii] +
                    f_flat[row0 + i + ii + 1] + f_flat[row1 + i + ii + 1]);
            }
        }

        i += 2;
    }

    // Trapezoidal rule for last column if nx is even
    if (nx % 2 == 0) {
        const double hx = x[nx - 1] - x[nx - 2];
        for (size_t j = 0; j + 1 < ny; ++j) {
            const double hy = y[j + 1] - y[j];
            const size_t row0 = j * nx;
            const size_t row1 = (j+1) * nx;
            total += 0.25 * hx * hy * (
                f_flat[row0 + nx - 2] + f_flat[row1 + nx - 2] +
                f_flat[row0 + nx - 1] + f_flat[row1 + nx - 1]);
        }
    }

    return total;
}

double simpson_2d_nonuniform_flat_weighted(
    const std::vector<double>& x,                        // full x vector, size nx
    const std::vector<double>& y,                        // full y vector, size ny
    const std::vector<double>& f_flat,                   // flattened f array, size nx * ny
    const std::vector<std::vector<double>>& Ax_weights,  // size (nx-2)/2 x 3
    const std::vector<std::vector<double>>& Ay_weights,  // size (ny-2)/2 x 3
    const std::vector<double>& dx,                       // size (nx-2)/2
    const std::vector<double>& dy                        // size (ny-2)/2
) {
    const size_t nx = x.size(); // pass this into integrator for slight speedup
    const size_t ny = y.size();

    if (f_flat.size() != nx * ny)
        throw std::invalid_argument("f_flat size must equal nx * ny");

    auto f = [&](size_t j, size_t i) -> double {
        return f_flat[j * nx + i];
    };

    double total = 0.0;

    // Bulk integration with precomputed Simpson weights
    for (size_t i = 0; i + 2 < nx; i += 2) {
        size_t wx_idx = i / 2;
        for (size_t j = 0; j + 2 < ny; j += 2) {
            size_t wy_idx = j / 2;

            double local = 0.0;
            for (int jj = 0; jj < 3; ++jj) {
                for (int ii = 0; ii < 3; ++ii) {
                    local += Ax_weights[wx_idx][ii] * Ay_weights[wy_idx][jj] * f(j + jj, i + ii);
                }
            }
            total += dx[wx_idx] * dy[wy_idx] * local / 4.0;
        }
    }

    // Handle last row (if ny even) with trapezoidal rule along y
    if (ny % 2 == 0) {
        size_t last_y = ny - 2;
        double hy = y[ny - 1] - y[last_y];

        for (size_t i = 0; i + 1 < nx; ++i) {
            double hx = x[i + 1] - x[i];

            total += 0.25 * hx * hy * (
                f(last_y, i) + f(last_y + 1, i) +
                f(last_y, i + 1) + f(last_y + 1, i + 1)
            );
        }
    }

    // Handle last column (if nx even) with trapezoidal rule along x
    if (nx % 2 == 0) {
        size_t last_x = nx - 2;
        double hx = x[nx - 1] - x[last_x];

        for (size_t j = 0; j + 1 < ny; ++j) {
            double hy = y[j + 1] - y[j];

            total += 0.25 * hx * hy * (
                f(j, last_x) + f(j + 1, last_x) +
                f(j, last_x + 1) + f(j + 1, last_x + 1)
            );
        }
    }

    return total;
}

double adaptive_simpson_recursive(const std::function<double(double)>& f,
                                  double a, double b,
                                  double fa, double fb, double fm,
                                  double eps, int depth, int max_depth) {
    double h = b - a;
    double c = (a + b) / 2.0;
    double fd = f((a + c) / 2.0);
    double fe = f((c + b) / 2.0);

    double Sleft  = (h / 12.0) * (fa + 4.0 * fd + fm);
    double Sright = (h / 12.0) * (fm + 4.0 * fe + fb);
    double S2 = Sleft + Sright;

    double S = (h / 6.0) * (fa + 4.0 * fm + fb); // coarse Simpson
    if (depth >= max_depth || std::abs(S2 - S) < 15.0 * eps) {
        return S2 + (S2 - S) / 15.0; // Richardson extrapolation
    } else {
        return adaptive_simpson_recursive(f, a, c, fa, fm, fd, eps / 2.0, depth + 1, max_depth)
             + adaptive_simpson_recursive(f, c, b, fm, fb, fe, eps / 2.0, depth + 1, max_depth);
    }
}

double adaptive_simpson(const std::function<double(double)>& f,
                        double a, double b, double eps,
                        int max_depth) {
    if (a == b) return 0.0;

    double fa = f(a);
    double fb = f(b);
    double fm = f((a + b) / 2.0);

    return adaptive_simpson_recursive(f, a, b, fa, fb, fm, eps, 0, max_depth);
}

double adaptive_simpson_2d(const std::function<double(double, double)>& f2d,
                           double x0, double x1,double y0, double y1,
                           double eps, int max_depth) {
    auto outer_integrand = [&](double x) {
        auto inner_integrand = [&](double y) {
            return f2d(x, y);
        };
        return adaptive_simpson(inner_integrand, y0, y1, eps / 2.0, max_depth);
    };

    return adaptive_simpson(outer_integrand, x0, x1, eps / 2.0, max_depth);
}

// Im(Si(x))=0 for real x
// Im(Ci(x))=pi for x<0 (0 otherwise)
double Si(double x) {
    if (x == 0.0) return 0.0;

    auto sin_integrand = [](double t) {
        return (t == 0.0) ? 1.0 : std::sin(t) / t;
    };

    const int nt = 200; // integration steps
    std::vector<double> t_vals = linspace(0.0, x, nt);
    std::vector<double> sin_integrand_vals(nt);
    
    for (int j = 0; j < nt; j++) {
        const auto t = t_vals[j];
        sin_integrand_vals[j] = sin_integrand(t);
    }

    return simpson_integrate(t_vals, sin_integrand_vals);
}

double Ci(double x) {
    if (x == 0.0) return -INFINITY;

    auto cos_integrand = [](double t) {
        if (std::abs(t) < 1e-4) return -0.5 * t; // Taylor exp near zero (otherwise simpson integration breaks)
        return (std::cos(t) - 1.0) / t;
    };

    const int nt = 200; // integration steps
    std::vector<double> t_vals = linspace(0.0, x, nt);
    std::vector<double> cos_integrand_vals(nt);
    
    for (int j = 0; j < nt; j++) {
        const auto t = t_vals[j];
        cos_integrand_vals[j] = cos_integrand(t);
    }

    return gamma_euler + std::log(abs(x)) + simpson_integrate(t_vals, cos_integrand_vals);
}

std::pair<double, double> SiCi(double x, const size_t n) {
    if (x == 0.0) return {0.0, -INFINITY};

    auto sin_integrand = [](double t) {
        return (t == 0.0) ? 1.0 : std::sin(t) / t;
    };

    auto cos_integrand = [](double t) {
        if (std::abs(t) < 1e-4) return -0.5 * t; // Taylor exp near zero (otherwise simpson integration breaks)
        return (std::cos(t) - 1.0) / t;
    };

    std::vector<double> t_vals = linspace(0.0, x, n);
    std::vector<double> sin_integrand_vals(n);
    std::vector<double> cos_integrand_vals(n);
    
    for (size_t j = 0; j < n; j++) {
        const auto t = t_vals[j];
        sin_integrand_vals[j] = sin_integrand(t);
        cos_integrand_vals[j] = cos_integrand(t);
    }

    const auto sin_int = simpson_integrate(t_vals, sin_integrand_vals);
    const auto cos_int = gamma_euler + std::log(abs(x)) + simpson_integrate(t_vals, cos_integrand_vals);

    return {sin_int, cos_int};
}

void read_sici_csv(const std::string& filename,
                   std::vector<double>& x_vals,
                   std::vector<double>& Si_vals,
                   std::vector<double>& Ci_vals)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    
    // Skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string x_str, si_str, ci_str;

        if (!std::getline(ss, x_str, ',')) continue;
        if (!std::getline(ss, si_str, ',')) continue;
        if (!std::getline(ss, ci_str, ',')) continue;

        try {
            double x = std::stod(x_str);
            double si = std::stod(si_str);
            double ci = std::stod(ci_str);

            x_vals.push_back(x);
            Si_vals.push_back(si);
            Ci_vals.push_back(ci);
        } catch (const std::exception& e) {
            std::cerr << "Warning: skipped invalid line: " << line << "\n";
        }
    }

    file.close();
    std::cout << "Read " << x_vals.size() << " rows from " << filename << "\n";
}

// std::pair<std::vector<double>, std::vector<double>> runge_kutta_ODE_solver(ODE_System system, const double t0, const double tf, const double y0, const size_t n) {
//     // write check for if tf < t0 -> do something to integrate backwards
//     const auto h = (tf - t0) / n; // step size

//     const auto N = syst.size();

//     // define solution vectors
//     std::vector<double> t_vals(n);
//     std::vector<double> y_vals(n);

//     // set initial conditions
//     t_vals[0] = t0;
//     y_vals[0] = y0;

//     for (int i = 0; i < n-1; i++) {
//         const auto t = t_vals[i]; // t_i
//         const auto y = y_vals[i]; // y_i

//         for (int j =0; j < N; j++) { // iterate over all ODEs
//             const auto k1 = dydx[j](t, y); // current step
//             const auto k2 = dydx[j](t + h / 2.0, y + h * k1 / 2.0);
//             const auto k3 = dydx[j](t + h / 2.0, y + h * k2 / 2.0);
//             const auto k4 = dydx[j](t + h, y + h * k3); // full step

//             y_vals[j][i+1] = y + (h / 6.0) * (k1 + 2.0 * (k2 + k3) + k4); // y_{i+1}
//         }
//         t_vals[i+1] = t + h; // t_{i+1}
//     }

//     return {t_vals, y_vals};
// }

// write my own/make this better
// make it suitable to pass in just dvdxi (rather than only {dvdxi, dwdxi})
std::pair<std::vector<double>, std::vector<state_type>> rk4_solver(
    const deriv_func& dydx,
    double x0,
    double xf,
    const state_type& y0,
    size_t n
) {
    assert(n >= 2 && "Number of steps must be at least 2.");
    // assert(abs(xf - x0) > 1e-10);

    const double h = (xf - x0) / static_cast<double>(n - 1);
    std::vector<double> x_vals(n);
    std::vector<state_type> y_vals(n);

    x_vals[0] = x0;
    y_vals[0] = y0;

    auto add = [](const state_type& a, const state_type& b) {
        state_type res(a.size());
        for (size_t i = 0; i < a.size(); ++i) res[i] = a[i] + b[i];
        return res;
    };

    auto scale = [](const state_type& v, double s) {
        state_type res(v.size());
        for (size_t i = 0; i < v.size(); ++i) res[i] = v[i] * s;
        return res;
    };

    double x = x0;
    state_type y = y0;

    for (size_t i = 1; i < n; ++i) {
        state_type k1 = dydx(x, y);
        state_type k2 = dydx(x + h / 2.0, add(y, scale(k1, h / 2.0)));
        state_type k3 = dydx(x + h / 2.0, add(y, scale(k2, h / 2.0)));
        state_type k4 = dydx(x + h, add(y, scale(k3, h)));

        state_type incr = add(
            add(scale(k1, 1.0), scale(k2, 2.0)),
            add(scale(k3, 2.0), scale(k4, 1.0))
        );
        y = add(y, scale(incr, h / 6.0));

        x += h;
        x_vals[i] = x;
        y_vals[i] = y;
    }

    return {x_vals, y_vals};
}

// modified bisection method only to be used for functions with exactly one root!
double root_finder(std::function<double(double)> f, double a, double b, double tol, int max_iter) {
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0.0) {
        throw std::runtime_error("Bisection method interval not bracketed!");
    }

    if (!std::isfinite(fa) || !std::isfinite(fb)) {
        throw std::runtime_error("f(a) or f(b) is not finite.");
    }

    for (int i = 0; i < max_iter; ++i) {
        double c = 0.5 * (a + b);
        double fc = f(c);

        if (!std::isfinite(fc)) {
            // throw std::runtime_error("f(c) is not finite during bisection.");
        }

        // Check convergence
        // since interval is not necessarily bracketed, need to remove second condition
        // otherwise it thinks it successfully found the root once it has gone through 
        // the entire interval
        if (std::abs(fc) < tol) {
            return c;
        }

        if (fa * fc < 0.0) {
        // if (fc > 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    // std::cout << "Bisection method failed, finding minimum of residual instead.\n";
    // return find_minimum(f, a, b);
    throw std::runtime_error("Bisection method did not converge.");
}

double root_finder2(std::function<double(double)> f, double a, double b, double tol, int max_iter) {
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0.0) {
        throw std::runtime_error("Bisection method interval not bracketed!");
    }

    if (!std::isfinite(fa) || !std::isfinite(fb)) {
        throw std::runtime_error("f(a) or f(b) is not finite.");
    }

    for (int i = 0; i < max_iter; ++i) {
        double c = 0.5 * (a + b);
        double fc = f(c);

        if (!std::isfinite(fc)) {
            // throw std::runtime_error("f(c) is not finite during bisection.");
        }

        // Check convergence
        // since interval is not necessarily bracketed, need to remove second condition
        // otherwise it thinks it successfully found the root once it has gone through 
        // the entire interval
        if (std::abs(fc) < tol) {
            return c;
        }

        if (fa * fc < 0.0) {
        // if (fc > 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    // std::cout << "Bisection method failed, finding minimum of residual instead.\n";
    // return find_minimum(f, a, b);
    throw std::runtime_error("Bisection method did not converge.");
}

std::vector<std::pair<double, double>> find_brackets(const std::function<double(double)>& f, double a, double b, int N) {
    std::vector<std::pair<double, double>> brackets;
    double dx = (b - a) / N;
    double x0 = a;
    double f0 = f(x0);

    for (int i = 1; i <= N; ++i) {
        double x1 = a + i * dx;
        double f1 = f(x1);
        if (f0 * f1 < 0) {
            brackets.emplace_back(x0, x1);
        }
        x0 = x1;
        f0 = f1;
    }

    return brackets;
}

std::vector<double> find_all_roots(
    const std::function<double(double)>& f,
    double a,
    double b,
    int N)
{
    auto brackets = find_brackets(f, a, b, N);
    std::vector<double> roots;
    for (const auto& [x0, x1] : brackets) {
        try {
            double root = root_finder(f, x0, x1);
            roots.push_back(root);
        } catch (...) {
            // skip if bisection fails
        }
    }
    return roots;
}

double find_smallest_root(
    const std::function<double(double)>& f,
    double a,
    double b,
    int N)
{
    auto roots = find_all_roots(f, a, b, N);
    if (roots.empty())
        throw std::runtime_error("No roots found.");
    return *std::min_element(roots.begin(), roots.end());
}

alglib::real_1d_array vector_to_real_1d_array(const std::vector<double>& vec) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i != vec.size() - 1) oss << ",";
    }
    oss << "]";
    
    alglib::real_1d_array result;
    result = oss.str().c_str();  // real_1d_array supports string assignment
    return result;
}