#include "deepphase.hpp"

int main(int argc, char** argv) {

    const PhaseTransition::Universe un;

    const PhaseTransition::PTParams params(0.5, 0.1, 1.0, 10.0, 1.71, "bag", "exp", un);

    const Hydrodynamics::FluidProfile profile(params);

    profile.write("DeepPhase/fluid_profile.csv");

    return 0;
}