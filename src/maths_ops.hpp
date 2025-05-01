// maths_ops.hpp
#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <initializer_list>

/*
TO DO:
- write logspace function (logarithmic version of linspace)
*/

template<typename T>
class vec {
public:
    // Constructors
    vec() = default;
    vec(size_t size) : data_(size, T{}) {}
    vec(std::initializer_list<T> list) : data_(list) {}

    // Size
    size_t size() const { return data_.size(); }

    // Access
    T& operator[](size_t i) { return data_[i]; }
    const T& operator[](size_t i) const { return data_[i]; }

    // Arithmetic
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

template <typename T>
class CubicSpline {
  public:
    CubicSpline();
    CubicSpline(const std::vector<T>& x, const std::vector<T>& y);

    T operator()(T xi) const;

  private:
    std::vector<T> x, y, h, a, b, c, d; // change these to have '_'?
    bool initialized_; // false if created by dflt ctor
};
#include "CubicSpline.tpp"

/**
 * @brief C++ version of numpy's linspace function
 *
 * @param start Initial value in the vector
 * @param end Final value in the vector
 * @param num Number of elements in the vector
 *
 * @return vec of length 'num' whose elements are linearly spaced between 'start' and 'end'
 */
std::vector<double> linspace(double start, double end, std::size_t num);

void print_vector(const std::vector<double> &vec);

double power(double x, int exp);

double power6(double x);