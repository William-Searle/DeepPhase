// PhaseTransition.cpp
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <unordered_set>
#include <cassert>

#include "phasetransition.hpp"

/*
TO DO:
- update default vals for beta, wp, wm in PTParams
- check if UF_trans is the same for vp and vm
- only inputs to PTParams should be model, alpha, vw and nucleation type i think
*/

namespace PhaseTransition {

// Universe
Universe::Universe()
    : Universe(dflt_universe::T0, dflt_universe::Ts, dflt_universe::H0, dflt_universe::Hs, dflt_universe::g0, dflt_universe::gs) {}

Universe::Universe(double T0, double Ts, double H0, double Hs, double g0, double gs)
    : T0_(T0),
      Ts_(Ts), 
      H0_(H0), 
      Hs_(Hs), 
      g0_(g0), 
      gs_(gs) {}

std::ostream& operator<<(std::ostream& os, const Universe& un) {
    os << "************** Universe parameters **************\n"
       << std::left
       << std::setw(20) << " " << std::setw(13) << "Today" << "Start of PT\n"
       << std::setw(20) << " " << std::setw(13) << "-----" << "-----------\n"
       << std::setw(20) << "Temperature:" << "T0=" << std::setw(10) << un.T0() << "Ts=" << un.Ts() << "\n"
       << std::setw(20) << "Hubble constant:" << "H0=" << std::setw(10) << un.H0() << "Hs=" << un.Hs() << "\n"
       << std::setw(20) << "Number of DoF:" << "g0=" << std::setw(10) << un.g0() << "gs=" << un.gs() << "\n"
       << "*************************************************\n";
       
    return os;
}

// some way to combine this with PTParams print()?
void Universe::print() const {
    std::cout << *this;
}

const Universe& default_universe() {
    static Universe u;;
    return u;
}

// PTParams
// make new ctor for reading in Veff since we calculate all the params
// alpha only exists for Bag model, maybe change how it is input/stored in PTParams
PTParams::PTParams()
    : PTParams(dflt_PTParams::vw, dflt_PTParams::alpha, dflt_PTParams::beta, dflt_PTParams::dtau, dflt_PTParams::wN, dflt_PTParams::model, dflt_PTParams::nuc_type, default_universe()) {}

PTParams::PTParams(double vw, double alpha, double beta, double dtau, double wN, const char* model, const char* nuc_type, const Universe& un)
    : vw_(vw),
      alpha_(alpha), 
      beta_(beta),
      dtau_(dtau),
      wN_(wN),
      tau_s_(),
      tau_fin_(),
      model_(),
      nuc_type_(),
      wall_type_(),
      Rs_(),
      vcj_(), cpsq_(), cmsq_(),
      universe_(un)
    {
      // defaults to bag model if input model is not in valid_models
      // write something that indicates other ctor should be called for Veff
      static const std::unordered_set<const char*> valid_models = {"bag"};//, "improved bag"};
      if (valid_models.count(model)) {
        model_ = model;
      } else {
        std::cout << "Warning: Invalid model '" << model << "' for equation of state. Using default model (" << dflt_PTParams::model << ")\n";
        model_ = dflt_PTParams::model;
      }

      // check valid bubble nucleation type
      static const std::unordered_set<const char*> valid_nuc = {"exp", "sim"};
      if (valid_nuc.count(nuc_type)) {
        nuc_type_ = nuc_type;
      } else {
        std::cout << "Warning: Invalid model '" << nuc_type << "' for bubble nucleation. Using default nucleation type (" << dflt_PTParams::nuc_type << ")\n";
        nuc_type_ = dflt_PTParams::nuc_type;
      }

      // check valid PT duration
      if (dtau < 0.0) {
        std::cout << "Warning: Invalid PT duration (dtau < 0). Using default value dtau=" << dflt_PTParams::dtau << "\n";
        dtau_ = dflt_PTParams::dtau;
      }
      tau_s_ = 1.0 / universe_.Hs();
      tau_fin_ = tau_s_ + dtau_;

      /***** calc model-dependent quantities *****/
      if (model == "bag") {
        cpsq_ = 1.0 / 3.0;
        cmsq_ = cpsq_;
      } else {
          throw std::invalid_argument("Only Bag model has been implemented so far");
      }

      Rs_ = std::pow(8 * M_PI, 1. / 3.) * vw_ / beta_;
      vcj_ = (1.0 / std::sqrt(3.0)) * (1.0 + std::sqrt(alpha_ + 3.0 * alpha_ * alpha_)) / (1.0 + alpha_);
      /*******************************************/

      /************** set wall type **************/
      const auto vcjsq = vcj_ * vcj_;
      const auto vwsq = vw_ * vw_;
      assert(vcjsq > cpsq_); // does this always hold?

      if (vwsq <= cpsq_) {
        wall_type_ = "deflagration";
      } else if (vwsq > cpsq_ && vwsq < vcjsq) {
        wall_type_ = "hybrid";
      } else {
        wall_type_ = "detonation";
      }
      /*******************************************/
    }

std::ostream& operator<<(std::ostream& os, const PTParams& params) {
    os << "********** Phase Transition parameters **********\n"
       << std::left
       << std::setw(35) << "Equation of state:" << params.model_ << "\n"
       << std::setw(35) << "Nucleation type:" << params.nuc_type_ << "\n"
       << std::setw(35) << "Wall type:" << params.wall_type_ << "\n"
       << std::setw(35) << "Wall velocity:" << "vw=" << params.vw_ << "\n"
       << std::setw(35) << "PT strength parameter:" << "alpha=" << params.alpha_ << "\n"
       << std::setw(35) << "Transition rate parameter:" << "beta=" << params.beta_ << "\n"
       << std::setw(35) << "PT duration" << "dtau=" << params.dtau_ << "\n"
       << std::setw(35) << "Speed of sound (broken phase):" << "cpsq=" << params.cpsq_ << "\n"
       << std::setw(35) << "Speed of sound (old phase):" << "cmsq=" << params.cmsq_ << "\n"
       << std::setw(35) << "Mean bubble separation:" << "Rs=" << params.Rs_ << "\n"
       << "*************************************************\n";
       
    return os;
}

void PTParams::print() const {
    std::cout << *this;
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

/********************/

} // namespace PhaseTransition