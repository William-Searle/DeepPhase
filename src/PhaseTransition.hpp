// PhaseTransition.hpp
#pragma once

#include <string>

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

        double get_cpsq() const { return cpsq_; } // speed of sound squared (symmetric phase)
        double get_cmsq() const { return cmsq_; } // speed of sound squared (broken phase)
        double get_vw() const { return vw_; } // wall velocity
        double get_alpha() const { return alpha_; } // strength parameter
        double get_beta() const { return beta_; } // inverse PT duration
        double get_Rs() const { return Rs_; } // characteristic length scale R_*
        std::string get_nuc_type() const { return nuc_type_; } // bubble nucleation type
      
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
        double get_T0() const { return T0_; } // temperature of universe
        double get_Ts() const { return Ts_; }

        double get_H0() const { return H0_; } // Hubble constant
        double get_Hs() const { return Hs_; }

        double get_g0() const { return g0_; } // number of dof
        double get_gs() const { return gs_; }
    
      private:
        const double T0_, Ts_, H0_, Hs_, g0_, gs_;
    };

} // namespace PhaseTransition