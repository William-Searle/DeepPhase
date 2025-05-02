// physics.cpp
#include <cmath>

#include "physics.hpp"

namespace Physics {

double gamma(double v) {
    return 1.0 / std::sqrt(1-std::pow(v,2));
}

double gammaSq(double v) {
    return 1.0 / (1-std::pow(v,2));
}

} // namespace Physics