// PhaseTransition.cpp
#include <string>
#include <cmath>
#include <iostream>
#include <vector>

#include "PhaseTransition.hpp"

/*
TO DO:
- update default vals for beta, wp, wm in PTParams
- check if UF_trans is the same for vp and vm
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

std::ostream& operator<<(std::ostream& os, const PTParams& params) {
    os << "cpsq = " << params.cpsq_ << "\n"
       << "cmsq = " << params.cmsq_ << "\n"
       << "vw = " << params.vw_ << "\n"
       << "alpha = " << params.alpha_ << "\n"
       << "beta = " << params.beta_ << "\n"
       << "nuc_type = " << params.nuc_type_ << "\n"
       << "Rs = " << params.Rs_ << "\n";
    return os;
}

// same for vp and vm?
double PTParams::vUF(const double v) const {
    return (vw_ - v) / ( 1.0 - vw_ * v);
}

/* Bag model */
// see (B.6) in Hindmarsh GWs from FOPT in the SSM
std::vector<double> PTParams::vpm() const { 
    std::vector<double> vpm;

    const auto vm = vw_;

    int sgn = vm > 1.0 / std::sqrt(3.0) ? 1 : -1;
    const auto tt = vm / 2.0 + 1.0 / (6.0 * vm);
    const auto troot = tt * tt + alpha_ * alpha_ + (2.0 / 3.0) * alpha_ - 1.0 / 3.0;
    const auto vp = (1.0 / (1.0 + alpha_)) * tt + sgn * std::sqrt(troot);

    vpm.push_back(vp);
    vpm.push_back(vm);

    return vpm; // {vp, vm}
}

std::vector<double> PTParams::wpm() const { 
    std::vector<double> wpm;

    const auto wp = 0.1; // PLACEHOLDER
    const auto wm = 0.1; // PLACEHOLDER

    wpm.push_back(wp);
    wpm.push_back(wm);

    return wpm; // {wp, wm}
}
/********************/

void PTParams::print() const {
    std::cout << *this;
}

// Universe
Universe::Universe()
    : Universe(2.725, 100.0, 67.8, 1.0, 3.91, 100.0) {}

Universe::Universe(double T0, double Ts, double H0, double Hs, double g0, double gs)
    : T0_(T0),
      Ts_(Ts), 
      H0_(H0), 
      Hs_(Hs), 
      g0_(g0), 
      gs_(gs) {}

} // namespace PhaseTransition