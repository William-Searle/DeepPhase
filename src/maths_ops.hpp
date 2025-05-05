// maths_ops.hpp

/**
 * @file maths_ops.hpp
 * @brief Mathematical utilities including a custom vector class and cubic spline interpolation.
 *
 * Provides a template `vec<T>` class with arithmetic operations, dot product, and norm.
 * Includes a `CubicSpline<T>` class for 1D spline interpolation and utility functions like `linspace`.
 */

#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <initializer_list>

/*
TO DO:
- write logspace function (logarithmic version of linspace)
- change cublicspline private vars to include '_' at end for consistency
*/

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

    T operator()(T xi) const;  // Evaluate spline at point xi

  private:
    std::vector<T> x_, y_;          // Input data
    std::vector<T> h_;              // Interval widths
    std::vector<T> a_, b_, c_, d_;  // Spline coefficients

    bool initialized_ = false;

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
std::vector<double> linspace(double start, double end, std::size_t num);

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

std::string to_string_with_precision(double value, int precision = 2);