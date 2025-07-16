<h1 align="center">
DeepPhase
</h1>

<div align="center">
<i>Some snarky catchphrase</i>
</div>

**DeepPhase** is a C++17 software package for calculating gravitational wave spectra from a first-order phase transition using either the Bag model or a via a direct calculation of thermodynamics from the effective potential.

## Dependencies

You need a C++17 compliant compiler, together with the following dependencies:
* CMake, version x.x.x or higher
* OpenMP, version x.x.x or higher
* GSL, version x.x or higher
* ALGLIB, version x.x or higher

On *Ubuntu/Debian*-based distributions, `ALGLIB` and `GSL` can be installed by running:

    sudo apt install libalglib-dev libgsl-dev

On *Fedora*-based distributions, instead use:

    sudo dnf install alglib-devel gsl-devel

Finally on *Mac*:

    brew install gsl alglib


Optionally, we support plotting using `matplotlib-cpp`. To enable this functionality, please install the `matplotlibcpp.h` header into `DeepPhase/includes` directory. This is then automatically detection during compilation and relevant plotting functionality is enabled.

## Building
To build the shared library and examples, use: *URL subject to change*

    git clone https://github.com/DeepPhase/DeepPhase
    cd DeepPhase
    mkdir build
    cd build
    cmake ..
    make

## Running
If the library and examples were successfully built, the examples and tests are available to run using the following executables: *subject to change*

    ./bin/run_fluid_profile
    ./bin/run_kinetic_spectrum
    ./bin/run_gw_spectrum
    ./bin/unit_tests

## Notes
* Bugs still present in GW calculation (not ready for use)
* Fluid profile & kinetic power spectrum working well!
