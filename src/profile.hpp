// profile.hpp
#pragma once

#include <vector>
#include <string>

#include "maths_ops.hpp"
#include "PhaseTransition.hpp"

/*
TO DO:
- pass in PTParams to FluidProfile and FluidSystem class and calculate initial state (y0) in ctor
- move generate_stream_plot() elsewhere? independent of initial conditions so the same for all instances of FluidProfile
- add some safety thing for private vales and prof in FluidProfile to stop them from changing (helper static function to precompute in initialiser list)
- make initial state (y0) const after calling FluidProfile ctor somehow?
- make copy of PTParams parameters in FluidProfile class
*/

namespace Hydrodynamics {

using state_type = std::vector<double>;
// using state_type = vec<double>;
using interp_type = CubicSpline<double>;

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

// not sure if this is the best way to implement fluid profile..
/**
 * @class FluidProfile
 * @brief Represents the hydrodynamic profile of a bubble wall in a first-order phase transition.
 */
class FluidProfile {
  public:
    FluidProfile(const PhaseTransition::PTParams& params); // ctor

    PhaseTransition::PTParams params() const { return params_; }; // PT parameters
    std::vector<double> init_state() const { return y0_; }; // Initial state vector {xi0, v0, w0}
    
    state_type xi_vals() const { return xi_vals_; }; // Vector of xi=r/t
    
    // Interpolation functions for profile
    state_type v_vals() const { return v_vals_; };
    interp_type v_prof() const { return v_prof_; };

    interp_type w_prof() const { return w_prof_; }; // not used, but might be good to keep for future
    state_type w_vals() const { return w_vals_; };

    interp_type la_prof() const { return la_prof_; };
    state_type la_vals() const { return la_vals_; };

    double v_shock(double xi) const; // v profile shock wave
    double w_shock(double xi) const; // w profile shock wave

    void write(const std::string& filename = "bubble_prof.csv") const; // write bubble profile to disk
    void plot(const std::string& filename = "bubble_prof.png") const; // Plots bubble profiles
    void generate_streamplot_data(int xi_pts=30, int y_pts=30, const std::string& filename="streamplot_data.csv") const;

  private:
    const PhaseTransition::PTParams params_;
    const double csq_;
    std::vector<double> y0_;
    state_type xi_vals_, v_vals_, w_vals_, la_vals_, la_vals_test_; // WARNING: Not const
    interp_type v_prof_, w_prof_, la_prof_, la_prof_test_;
    // NOTE: removed 'const' for vals and prof so they can be defined in the ctor body

    // put number of integration points in input file? seems bad to hardcode
    /**
     * @brief Solves the hydrodynamic profile ODE system.
     * 
     * @param n Number of integration steps.
     * 
     * @return Vector of state vectors representing the profile.
     */
    std::vector<state_type> solve_profile(int n=100) const;

    /**
     * @brief Constructs cubic splines for the hydrodynamic profiles.
     * 
     * @return Vector of cubic splines {v, w, lambda}.
     */
    void profile(bool read_prof = true);

    std::vector<state_type> read(const std::string& filename) const; // read bubble profile from disk

    state_type calc_lambda_vals(state_type w_vals) const;
    interp_type calc_lambda_prof() const;
};

} // namespace Hydrodynamics