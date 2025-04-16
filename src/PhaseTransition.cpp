// PhaseTransition.cpp
#include <string>
#include <cmath>

#include "constants.hpp"
#include "PhaseTransition.hpp"

/*
TO DO:
- update default val for beta in PTParams
- better way to define Rs in PTParams? (currently have to write definition twice - one for each constructor)
*/

namespace PhaseTransition {

// PTParams
PTParams::PTParams()
    : cpsq_(1. / 3.), cmsq_(1. / 3.), vw_(0.7), alpha_(0.1), beta_(1.0),
      Rs_(std::pow(8 * PI, 1. / 3.) * vw_ / beta_), nuc_type_("exp") {}

PTParams::PTParams(double cpsq, double cmsq, double vw, double alpha, double beta, std::string nuc_type)
    : cpsq_(cpsq), cmsq_(cmsq), vw_(vw), alpha_(alpha), beta_(beta), 
      Rs_(std::pow(8*PI, 1./3.) * vw_ / beta_), nuc_type_(nuc_type) {}

// Universe
Universe::Universe()
    : T0_(2.725), Ts_(100.0), H0_(67.8), Hs_(1.0), g0_(3.91), gs_(100.0) {}

Universe::Universe(double T0, double Ts, double H0, double Hs, double g0, double gs)
    : T0_(T0), Ts_(Ts), H0_(H0), Hs_(Hs), g0_(g0), gs_(gs) {}

} // namespace PhaseTransition