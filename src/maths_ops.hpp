// maths_ops.hpp
#pragma once

/*
TO DO:
- write logspace function (logarithmic version of linspace)
*/

#include <vector>

/**
 * @brief C++ version of numpy's linspace function
 *
 * @param start Initial value in the vector
 * @param end Final value in the vector
 * @param num Number of elements in the vector
 *
 * @return Vector of length 'num' whose elements are linearly spaced between 'start' and 'end'
 */
std::vector<double> linspace(double start, double end, std::size_t num);