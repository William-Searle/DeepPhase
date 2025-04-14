// main.cpp
#include <iostream>
#include <array>

#include "constants.hpp"
#include "PhaseTransition.hpp"
#include "hydrodynamics.hpp"
#include "spectrum.hpp"

/*
TO DO:
-separate into different test files to test each function so I can compare to the python output
-use libraries?
*/

int main() {
    /* hydrodynamics test */
    // const auto xi = 0.5;
    // const auto v = 0.6;
    // const auto w = 0.5;
    // const auto T = 100.0;
    // const auto cssq = 1./3.;
    // const std::array<double,3> xiw = {xi,w,T};

    // const auto df_dv = Hydrodynamics::dfdv(xiw, v, cssq);

    // std::cout << "mu: " << Hydrodynamics::mu(1./2., 1./3.) << "\n"
    //           << "getwow: " << Hydrodynamics::getwow(1./2., 1./3.) << "\n"
    //           << "dfdv: ";
    // for (const auto &val : df_dv) {
    //     std::cout << val << " ";
    // }
    // std::cout << "\n";

    // PhaseTransition classes test
    const PhaseTransition::Universe u;
    const PhaseTransition::PTParams params;
    std::cout << "(T0,Ts,H0,Hs,g0,gs)=(" << u.get_T0() << ","
                                         << u.get_Ts() << ","
                                         << u.get_H0() << ","
                                         << u.get_Hs() << ","
                                         << u.get_g0() << ","
                                         << u.get_gs() << ")"
                                         << "\n";
    std::cout << "(cpsq,cmsq,vw,alpha)=(" << params.get_cpsq() << ","
                                          << params.get_cmsq() << ","
                                          << params.get_vw() << ","
                                          << params.get_alpha() << ")"
                                          << "\n";

    return 0;
}
