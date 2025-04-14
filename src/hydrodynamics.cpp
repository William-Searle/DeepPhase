// hydrodynamics.cpp

/*
TO DO:
- change short functions to inline functions and move to header
*/

#include <cmath>
#include <array>

namespace Hydrodynamics {

double mu(double xi, double v) {
    return (xi - v) / (1.0 - xi * v);
}

double getwow(double vp, double vm) {
    const auto wp = vp / (1. - std::pow(vp, 2));
    const auto wm = vm / (1. - std::pow(vm, 2));

    return wp/wm;
}

std::array<double,3> dfdv(const std::array<double,3>& xiw, double v, double csq) {
    if (xiw.size() != 3) {
        throw std::invalid_argument("xiw must be an array of the form {xi, w, T}");
    }
    const auto& [xi, w, T] = xiw;
    const auto dxidv  = (std::pow(mu(xi, v),2) / csq - 1) * (1 - v*xi) 
                        * (xi/2) * (1/v) / (1-std::pow(v,2));             // velocity profile
    const auto dwdv = (1 + 1/csq) * mu(xi, v) * w / (1 - std::pow(v,2)); // enthalpy profile
    const auto dTdv = mu(xi,v) * T / (1 - std::pow(v,2));                 // temperature profile

    // change this and xiw to vec?
    // set up template function so output matches data type of xiw?
    const std::array<double,3> derivs = {dxidv, dwdv, dTdv};

    return derivs;
}
}