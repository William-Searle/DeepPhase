// maths_ops.cpp
#include <vector>
#include <iostream>
#include <iomanip>

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