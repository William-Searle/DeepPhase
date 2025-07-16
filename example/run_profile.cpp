#include "deepphase.hpp"

int main() {

    const PhaseTransition::Universe un;

    const PhaseTransition::PTParams params1(0.5, 0.1, 1.0, 10.0, 1.71, "bag", "exp", un);
    const Hydrodynamics::FluidProfile profile1(params1);

    const auto kRs_vals = logspace(1e-1, 1e+3, 500);

    const auto Ek1 = Spectrum::Ekin(kRs_vals, profile1);
    const auto Eks1 = Spectrum::zetaKin(Ek1);

    std::cout << "test run" << std::endl;
}