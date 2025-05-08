// PhaseTransition.cpp
#include <string>
#include <cmath>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <cassert>

#include "PhaseTransition.hpp"

/*
TO DO:
- update default vals for beta, wp, wm in PTParams
- check if UF_trans is the same for vp and vm
- only inputs to PTParams should be model, alpha, vw and nucleation type i think
*/

namespace PhaseTransition {

// PTParams
PTParams::PTParams()
    : PTParams(0.5, 0.1, 1.0, "bag", "exp") {}

PTParams::PTParams(double vw, double alpha, double beta, std::string model, std::string nuc_type)
    : vw_(vw),
      alpha_(alpha), 
      beta_(beta),
      model_(model),
      nuc_type_(nuc_type),
      Rs_(std::pow(8 * M_PI, 1. / 3.) * vw_ / beta_),
      cpsq_(), cmsq_(),
      wall_type_()
    {
      // defaults to bag model if input model is not in valid_models
      static const std::unordered_set<std::string> valid_models = {"bag", "improved bag", "Veff"};
      model_ = valid_models.count(model) ? model : "bag";

      // compute model-dependent quantities
      const auto mod_par = model_params(model_);
      cpsq_ = mod_par[0];
      cmsq_ = mod_par[1];
      //   alpha_ = params[2];

      vcj_ = (1.0 / std::sqrt(3.0)) * (1.0 + std::sqrt(alpha_ + 3.0 * alpha_ * alpha_)) / (1.0 + alpha_);
      const auto vcjsq = vcj_ * vcj_;
      assert(vcjsq > cpsq_); // does this always hold?

      const auto vwsq = vw_ * vw_;

      if (vwsq <= cpsq_) {
        wall_type_ = "deflagration";
      } else if (vwsq > cpsq_ && vwsq < vcjsq) {
        wall_type_ = "hybrid";
      } else {
        wall_type_ = "detonation";
      }
    }

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

std::vector<double> PTParams::model_params(const std::string& model) const {
    if (model == "bag") {
        return bag_params();
    } else if (model == "improved bag") {
        throw std::invalid_argument("Only Bag model has been implemented so far");
        // return improved_bag_params();
    } else if (model == "Veff") {
        throw std::invalid_argument("Only Bag model has been implemented so far");
        // return veff_params();
    } else {
        std::cout << "Warning: Invalid equation of state model. Defaulting to Bag model." << std::endl;
        return bag_params();
    }
}

double PTParams::vUF(const double v) const { // universe frame (centre of bubble)
    return (vw_ - v) / ( 1.0 - vw_ * v);
}

/* Bag model: */
/*
p+ = (1/3) * a+ * T+^4 - eps
e+ = a+ * T+^4 + eps
p- = (1/3) * a- * T-^4
e- = a- * T-^4

csq = dp/de -> cpsq = cmsq = 1/3
alpha = eps / (a+ * T+^4)
*/

std::vector<double> PTParams::bag_params() const {
    std::vector<double> param_vec;
    // const auto eps =  // Bag constant
    const auto cpsq = 1.0 / 3.0;
    const auto cmsq = cpsq;
    // const auto alpha = 

    param_vec.push_back(cpsq);
    param_vec.push_back(cmsq);
    // param_vec.push_back(alpha);
    return param_vec; // {cpsq, cmsq, alpha}
}

// see (B.6) in Hindmarsh GWs from FOPT in the SSM
// only true for bag model - move all bag stuff together at some point
std::vector<double> PTParams::vpm() const {  // rest frame of bubble wall
    std::vector<double> vpm;

    const auto vm = vw_; // not sure if correct - see B.7

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