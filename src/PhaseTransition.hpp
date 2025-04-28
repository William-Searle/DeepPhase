// PhaseTransition.hpp
#pragma once

#include <string>
#include <iostream>

/*
TO DO:
- add option to input a Veff -> derive fluid dynamics from this?
- change PTParams to single ctor with default arguments (see slide 23, wk 4)
*/

namespace PhaseTransition {

    class PTParams { // will probably need to update this later
      public:
        // ctors
        PTParams();
        PTParams(double cpsq, double cmsq, double vw, double alpha, double beta, std::string nuc_type);

        double cpsq() const { return cpsq_; } // speed of sound squared (symmetric phase)
        double cmsq() const { return cmsq_; } // speed of sound squared (broken phase)
        double vw() const { return vw_; } // wall velocity
        double alpha() const { return alpha_; } // strength parameter
        double beta() const { return beta_; } // inverse PT duration
        double Rs() const { return Rs_; } // characteristic length scale R_*
        std::string nuc_type() const { return nuc_type_; } // bubble nucleation type

        // print params
        void print() const;
        friend std::ostream& operator<<(std::ostream& os, const PTParams& p);
      
      private:
          const double cpsq_, cmsq_, vw_, alpha_, beta_, Rs_;
          const std::string nuc_type_;
    };

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