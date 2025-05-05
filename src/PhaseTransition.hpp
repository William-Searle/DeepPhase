// PhaseTransition.hpp
#pragma once

#include <string>
#include <iostream>
#include <vector>

/*
TO DO:
- add option to input a Veff -> derive fluid dynamics from this?
- change PTParams to single ctor with default arguments (see slide 23, wk 4)
- initialise all the new variables i've added in ctor
- add deflag, hybrid, detonation type in PTParams
- add vp, vm, wp, wm as parameters in PTParams? currently uses function to calculate since they depend on Veff
  - these are calculated from a specific matching condition using bag model - generalise to other models in future
*/

namespace PhaseTransition {

/**
 * @class PTParams
 * @brief Class representing the parameters of the phase transition (PT).
 * 
 * This class encapsulates the parameters needed to describe the phase transition dynamics,
 * including speeds of sound, wall velocity, strength of the transition, and bubble nucleation type.
 */
class PTParams { // will probably need to update this later
  public:
    // ctors
    PTParams();
    PTParams(double cpsq, double cmsq, double vw, double alpha, double beta, std::string nuc_type);

    double cpsq() const { return cpsq_; } // speed of sound squared (symmetric phase)
    double cmsq() const { return cmsq_; } // speed of sound squared (broken phase)
    double csq() const { return cmsq_; } // keep this?
    double vw() const { return vw_; } // wall velocity
    double alpha() const { return alpha_; } // strength parameter
    double beta() const { return beta_; } // inverse PT duration
    double Rs() const { return Rs_; } // characteristic length scale R_*

    std::string nuc_type() const { return nuc_type_; } // bubble nucleation type

    double vUF(const double v) const; // converts v,w,e to universe frame

    std::vector<double> vpm() const; // fluid velocity {vp=symmetric phase, vm=broken phase}
    std::vector<double> wpm() const; // fluid enthalpy {wp=symmetric phase, wm=broken phase}

    // print params
    void print() const;
    friend std::ostream& operator<<(std::ostream& os, const PTParams& p);
  
  private:
      const double cpsq_, cmsq_, vw_, alpha_, beta_, Rs_;
      const std::string nuc_type_;
};

/**
 * @class Universe
 * @brief Class representing the universe parameters used in the phase transition calculations.
 * 
 * This class holds the current temperature, Hubble constant, and degrees of freedom of the universe.
 */
class Universe {
  public:
    // ctors
    Universe();
    Universe(double T0, double Ts, double H0, double Hs, double g0, double gs);
    
    // params today (0) and at start of PT (s)
    double T0() const { return T0_; } // temperature of universe
    double Ts() const { return Ts_; }

    double H0() const { return H0_; } // Hubble constant
    double Hs() const { return Hs_; }

    double g0() const { return g0_; } // number of dof
    double gs() const { return gs_; }

    void print() const; // print universe params

  private:
    const double T0_, Ts_, H0_, Hs_, g0_, gs_;
};

} // namespace PhaseTransition