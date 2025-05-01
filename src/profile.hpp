#pragma once

#include <vector>
#include <string>
#include "maths_ops.hpp"

/*
TO DO:
- pass in PTParams to FluidProfile and FluidSystem class and calculate initial state (y0) in ctor
- move generate_stream_plot() elsewhere? independent of initial conditions so the same for all instances of FluidProfile
*/

namespace Hydrodynamics {

using state_type = std::vector<double>;

class FluidSystem {
  public:
    FluidSystem(double csq) : csq_(csq) {}

    double dxi_dtau(double xi, double v) const;
    double dv_dtau(double xi, double v) const;
    double dw_dtau(double xi, double v, double w) const;

    void operator()(const state_type &y, state_type &dydt, double /*tau*/) const;

  private:
    double csq_;
};

struct push_back_state {
    std::vector<state_type> &states;
    std::vector<double> &times;

    push_back_state(std::vector<state_type> &s, std::vector<double> &t) : states(s), times(t) {}

    void operator()(const state_type& y, double t) const;
};

// not sure if this is the best way to implement fluid profile..
class FluidProfile {
  public:
    FluidProfile(state_type& y0, double csq); // ctor

    state_type init_state() const { return y0_; };
    double csq() const { return csq_; };
    std::vector<double> xi_vals() const { return xi_vals_; };
    CubicSpline<double> v_prof() const { return v_prof_; };
    CubicSpline<double> w_prof() const { return w_prof_; }; // not used, but might be good to keep for future
    CubicSpline<double> la_prof() const { return la_prof_; };

    void write() const; // not finished
    void generate_streamplot_data(int xi_pts=30, int y_pts=30, const std::string& filename="streamplot_data.csv") const;

  private:
    const double csq_;
    const state_type y0_; // initial state vector {xi0, v0, w0}
    std::vector<double> xi_vals_, v_vals_, w_vals_, la_vals_;
    CubicSpline<double> v_prof_, w_prof_, la_prof_;
    // NOTE: removed 'const' for vals and prof so they can be defined in the ctor body
    // add some safety thing to stop them from changing (helper static function to precompute in initialiser list)

    // put number of integration points in input file? seems bad to hardcode
    std::vector<state_type> solve_prof(int n=100) const;
    std::vector<CubicSpline<double>> profile();
};

} // namespace Hydrodynamics