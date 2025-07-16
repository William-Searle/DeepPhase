// maths_ops.hpp
#ifndef INCLUDE_MATHS_OPS_HPP_H
#define INCLUDE_MATHS_OPS_HPP_H

/**
 * @file maths_ops.hpp
 * @brief Mathematical utilities including a custom vector class and cubic spline interpolation.
 *
 * Provides a template `vec<T>` class with arithmetic operations, dot product, and norm.
 * Includes a `CubicSpline<T>` class for 1D spline interpolation and utility functions like `linspace`.
 */


#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <initializer_list>
#include <functional>

#include "ap.h"

// get rid of this
/**
 * @brief A custom vector class supporting arithmetic operations and linear algebra utilities.
 * 
 * @tparam T Type of the vector elements (typically double).
 */
template<typename T>
class vec {
public:
    // Constructors
    vec() = default;
    vec(size_t size) : data_(size, T{}) {}
    vec(std::initializer_list<T> list) : data_(list) {}

    // Returns size of vector
    size_t size() const { return data_.size(); }

    // Accesses the i-th element (modifiable & read-only).
    T& operator[](size_t i) { return data_[i]; }
    const T& operator[](size_t i) const { return data_[i]; }

    // Vector Arithmetic
    vec<T> operator+(const vec<T>& other) const;
    vec<T> operator-(const vec<T>& other) const;
    vec<T> operator*(T scalar) const;
    vec<T> operator/(T scalar) const;

    vec<T>& operator+=(const vec<T>& other);
    vec<T>& operator-=(const vec<T>& other);
    vec<T>& operator*=(T scalar);
    vec<T>& operator/=(T scalar);

    // Dot product
    T dot(const vec<T>& other) const;

    // Norm
    T norm() const;

    // Output
    void print() const;

private:
    std::vector<T> data_;
};
#include "vector.tpp"

/**
 * @brief Cubic spline interpolation class.
 * 
 * @tparam T Type for input/output values (typically double).
 */
template <typename T>
class CubicSpline {
  static_assert(std::is_floating_point<T>::value, "CubicSpline only supports floating-point types.");

  public:
    CubicSpline();  // Default constructor
    CubicSpline(const std::vector<T>& x, const std::vector<T>& y);  // Construct and compute spline
    CubicSpline(const vec<T>& x, const vec<T>& y); // for custom vector type

    void build(const std::vector<T>& x, const std::vector<T>& y); // Build splines
    bool is_initialised() const {return initialised_; } // check if spline is initialised
    void check_convergence() const; // check convergence

    T operator()(T xi) const;  // Evaluate spline at point xi

    // Scalar operations from the right
    CubicSpline<T> operator+(T scalar) const;
    CubicSpline<T> operator-(T scalar) const;
    CubicSpline<T> operator*(T scalar) const;
    CubicSpline<T> operator/(T scalar) const;

    // Scalar operations from the left
    template <typename U>
    friend CubicSpline<U> operator+(U scalar, const CubicSpline<U>& spline);

    template <typename U>
    friend CubicSpline<U> operator-(U scalar, const CubicSpline<U>& spline);

    template <typename U>
    friend CubicSpline<U> operator*(U scalar, const CubicSpline<U>& spline);

  private:
    std::vector<T> x_, y_;          // Input data
    std::vector<T> h_;              // Interval widths
    std::vector<T> a_, b_, c_, d_;  // Spline coefficients

    const T tol_ = static_cast<T>(1e-4); // convergence tolerance
    bool initialised_ = false;

    bool is_strictly_monotonic(const std::vector<T>& x) const; // Check monotonicity
};
#include "CubicSpline.tpp"

/**
 * @brief Creates a linearly spaced vector between `start` and `end` (C++ equivalent of python's linspace function).
 *
 * @param start Initial value.
 * @param end Final value.
 * @param num Number of points.
 * 
 * @return std::vector<double> Linearly spaced values.
 */
std::vector<double> linspace(double start, double end, std::size_t num=100);

std::vector<double> logspace(double start, double stop, std::size_t num=100);

/**
 * @brief Computes x raised to an integer exponent.
 * 
 * @param x Base value.
 * @param exp Integer exponent.
 * 
 * @return double Result of x^exp.
 * 
 */
double power(double x, int exp);

/**
 * @brief Computes x^6
 * 
 * @param x Base value.
 * 
 * @return double Result of x^6.
 */
double power6(double x);

double power3(double x);

std::string to_string_with_precision(double value, int precision = 2);

struct SimpsonWeights2D {
    std::vector<std::vector<double>> Ax_weights; // size: (nx-2) x 3
    std::vector<std::vector<double>> Ay_weights; // size: (ny-2) x 3
    std::vector<double> dx;  // size: nx-2
    std::vector<double> dy;  // size: ny-2
};

// static void precompute_1d_weights(
//     const std::vector<double>& coords,
//     std::vector<std::vector<double>>& weights,
//     std::vector<double>& intervals);

SimpsonWeights2D precompute_simpson_weights_2d(
    const std::vector<double>& x,
    const std::vector<double>& y);

double simpson_2d_nonuniform_flat_weighted(
    const std::vector<double>& x,                        // full x vector, size nx
    const std::vector<double>& y,                        // full y vector, size ny
    const std::vector<double>& f_flat,                   // flattened f array, size nx * ny
    const std::vector<std::vector<double>>& Ax_weights,  // size (nx-2)/2 x 3
    const std::vector<std::vector<double>>& Ay_weights,  // size (ny-2)/2 x 3
    const std::vector<double>& dx,                       // size (nx-2)/2
    const std::vector<double>& dy                        // size (ny-2)/2
);

double simpson_integrate(const std::vector<double>& x, const std::vector<double>& y);
double simpson_nonuniform(const std::vector<double>& x, const std::vector<double>& y);
double simpson_2d_integrate(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f);
double simpson_2d_integrate_flat(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& f_flat);
double simpson_2d_nonuniform_flat(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& f_flat);

double adaptive_simpson_recursive(const std::function<double(double)>& f,
                                  double a, double b,
                                  double fa, double fb, double fm,
                                  double eps, int depth, int max_depth);

double adaptive_simpson(const std::function<double(double)>& f,
                        double a, double b,
                        double eps = 1e-8, int max_depth = 20);

double adaptive_simpson_2d(const std::function<double(double, double)>& f2d,
                           double x0, double x1, double y0, double y1,
                           double eps = 1e-8, int max_depth = 20);

double Si(double x);
double Ci(double x);
std::pair<double, double> SiCi(double x, const size_t n=1000);

void read_sici_csv(const std::string& filename,
                   std::vector<double>& x_vals,
                   std::vector<double>& Si_vals,
                   std::vector<double>& Ci_vals);

using state_type = std::vector<double>;
using deriv_func = std::function<state_type(double, const state_type&)>;

std::pair<std::vector<double>, std::vector<state_type>> rk4_solver(
    const deriv_func& dydx,
    double x0,
    double xf,
    const state_type& y0,
    size_t n
);

double root_finder(std::function<double(double)> f, double a, double b, double tol = 1e-8, int max_iter = 100);

std::vector<std::pair<double, double>> find_brackets(const std::function<double(double)>& f, double a, double b, int N = 1000);
std::vector<double> find_all_roots(
    const std::function<double(double)>& f,
    double a,
    double b,
    int N = 1000);
double find_smallest_root(
    const std::function<double(double)>& f,
    double a,
    double b,
    int N = 1000);

alglib::real_1d_array vector_to_real_1d_array(const std::vector<double>& vec);

#endif // INCLUDE_MATHS_OPS_HPP_H
