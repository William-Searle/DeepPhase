// PhaseTransition.cpp
#include "PhaseTransition.hpp"

namespace PhaseTransition {

// PTParams
PTParams::PTParams()
    : cpsq_(1./3.), cmsq_(1./3.), vw_(0.7), alpha_(0.1) {}

PTParams::PTParams(double cpsq, double cmsq, double vw, double alpha)
    : cpsq_(cpsq), cmsq_(cmsq), vw_(vw), alpha_(alpha) {}

// Universe
Universe::Universe()
    : T0_(2.725), Ts_(100.0), H0_(67.8), Hs_(1.0), g0_(3.91), gs_(100.0) {}

Universe::Universe(double T0, double Ts, double H0, double Hs, double g0, double gs)
    : T0_(T0), Ts_(Ts), H0_(H0), Hs_(Hs), g0_(g0), gs_(gs) {}

} // namespace PhaseTransition