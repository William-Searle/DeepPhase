// profile.hpp
#pragma once

#include <vector>
#include <string>

#include "maths_ops.hpp"
#include "PhaseTransition.hpp"

/*
TO DO:
- move generate_stream_plot() elsewhere? independent of initial conditions so the same for all instances of FluidProfile
- add some safety thing for private vals and prof in FluidProfile to stop them from changing (helper static function to precompute in initialiser list)
- add dflt ctor for FluidProfile that just uses dflt ctor for PTParams
*/

namespace Hydrodynamics {

/**
 * @brief Computes the Lorentz factor between the wall frame and the universe frame.
 *
 * @param xi Fluid velocity in the wall frame.
 * @param v Fluid velocity in the universe frame.
 * 
 * @return Lorentz factor (mu).
 */
double mu(double xi, double v);

double dvdxi(double xi, double v, const double csq);
double dwdxi(double xi, double v, double w, const double csq);

using state_type = std::vector<double>;

/**
 * @class FluidSystem
 * @brief Defines the ODE system governing the evolution of fluid variables.
 */
class FluidSystem {
  public:
    FluidSystem(const PhaseTransition::PTParams& params) : params_(params) {}

    double dxi_dtau(double xi, double v) const;
    double dv_dtau(double xi, double v) const;
    double dw_dtau(double xi, double v, double w) const;

     /**
     * @brief Operator overload to evaluate the full ODE system at a given state and time.
     * 
     * @param y Current state vector {xi, v, w}.
     * @param dydt Output derivative vector.
     * @param tau Current value of the independent variable (not used explicitly).
     */
    void operator()(const state_type &y, state_type &dydt, double /*tau*/) const;

  private:
    PhaseTransition::PTParams params_;
};

/// @brief Observer to store integration states and times during ODE solving.
struct push_back_state {
    std::vector<state_type> &states;
    std::vector<double> &times;

    push_back_state(std::vector<state_type> &s, std::vector<double> &t) : states(s), times(t) {}

    void operator()(const state_type& y, double t) const;
};

/**
 * @class FluidProfile
 * @brief Represents the hydrodynamic profile of a bubble wall in a first-order phase transition.
 */
class FluidProfile {
  public:
    FluidProfile(const PhaseTransition::PTParams& params, const bool test=false); // ctor
    // write dflt ctor?

    // are these needed?
    PhaseTransition::PTParams params() const { return params_; }; // PT parameters
    std::vector<double> init_state() const { return y0_; }; // Initial state vector {xi0, v0, w0}
    
    state_type xi_vals() const { return xi_vals_; }; // Vector of xi=r/t
    state_type v_vals() const { return v_vals_; }; // v(xi)
    state_type w_vals() const { return w_vals_; }; // w(xi)
    state_type la_vals() const { return la_vals_; }; // la(xi)

    double v_shock(double xi) const; // v profile shock wave
    double w_shock(double xi) const; // w profile shock wave

    void write(const std::string& filename = "bubble_prof.csv") const; // write bubble profile to disk
    void plot(const std::string& filename = "bubble_prof.png") const; // Plots bubble profiles
    void generate_streamplot_data(int xi_pts=30, int y_pts=30, const std::string& filename="streamplot_data.csv") const;

  private:
    const PhaseTransition::PTParams params_; // local copy of PT parameters
    const double csq_; // local copy of speed of sound
    double xi0_, xif_;
    std::vector<double> y0_; // initial conditions
    state_type xi_vals_, v_vals_, w_vals_, la_vals_; // xi, v(xi), w(xi), la(x)
    int mode_;
    double vp_, vpUF_, vm_, vmUF_; // v+, v- in bubble wall frame and universe frame (UF)
    double wpwN_, wmwN_; // w+/wN, w-/wN
    double alpha_;

    void get_mode_y0(); // define hydrodynamic mode and initial state (xi0,v0,w0)
    double calc_vm(double vp) const;
    double calc_vp(double vm) const;
    double calc_wm(double wp, double vp, double vm) const;
    // put number of integration points in input file? seems bad to hardcode
    /**
     * @brief Solves the hydrodynamic profile ODE system.
     * 
     * @param n Number of integration steps.
     * 
     * @return Vector of state vectors representing the profile.
     */
    std::vector<state_type> solve_profile(int n=100) const;
    std::vector<state_type> read(const std::string& filename) const; // read bubble profile from disk
    state_type calc_lambda_vals(state_type w_vals) const;
};

} // namespace Hydrodynamics