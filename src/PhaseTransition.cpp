// PhaseTransition.cpp
#include <string>
#include <cmath>
#include <iostream>

#include "constants.hpp"
#include "PhaseTransition.hpp"

/*
TO DO:
- update default val for beta in PTParams
*/

namespace PhaseTransition {

// PTParams
PTParams::PTParams()
    : PTParams(1./3., 1./3., 0.7, 0.1, 1.0, "exp") {}

PTParams::PTParams(double cpsq, double cmsq, double vw, double alpha, double beta, std::string nuc_type)
    : cpsq_(cpsq), 
      cmsq_(cmsq), 
      vw_(vw), 
      alpha_(alpha), 
      beta_(beta), 
      nuc_type_(nuc_type), 
      Rs_(std::pow(8 * M_PI, 1. / 3.) * vw_ / beta_) {}

std::ostream& operator<<(std::ostream& os, const PTParams& p) {
    os << "cpsq = " << p.cpsq_ << "\n"
       << "cmsq = " << p.cmsq_ << "\n"
       << "vw = " << p.vw_ << "\n"
       << "alpha = " << p.alpha_ << "\n"
       << "beta = " << p.beta_ << "\n"
       << "nuc_type = " << p.nuc_type_ << "\n"
       << "Rs = " << p.Rs_ << "\n";
    return os;
}

void PTParams::print() const {
    std::cout << *this;
}

// Universe
// write dflt ctor same way as PTParam? Seems easier to read this way
Universe::Universe()
    : T0_(2.725), Ts_(100.0), H0_(67.8), Hs_(1.0), g0_(3.91), gs_(100.0) {}

Universe::Universe(double T0, double Ts, double H0, double Hs, double g0, double gs)
    : T0_(T0),
      Ts_(Ts), 
      H0_(H0), 
      Hs_(Hs), 
      g0_(g0), 
      gs_(gs) {}

} // namespace PhaseTransition