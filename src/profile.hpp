#pragma once

/*
TO DO:
- create class for fluid profile to  write to disk, streamplot, print (maybe just have this in vector class?)
*/

namespace FluidProfile {

using state_type = std::vector<double>;

class FluidSystem {
  public:
    FluidSystem(double csq) : csq_(csq) {}

    void operator()(const state_type &y, state_type &dydt, double /*tau*/) const;

  private:
    double csq_;

    double dxi_dtau(double xi, double v, double csq) const;
    double dv_dtau(double xi, double v, double csq) const;
    double dw_dtau(double xi, double v, double w, double csq) const;
};

struct push_back_state {
    std::vector<state_type> &states;
    std::vector<double> &times;

    push_back_state(std::vector<state_type> &s, std::vector<double> &t) : states(s), times(t) {}

    void operator()(const state_type &y, double t) const;
};

void solve_prof(state_type y0, double csq, int n=100);
void plot_velocity_profile(const std::string& filename = "fluid_solution.csv");

} // namespace FluidProfile