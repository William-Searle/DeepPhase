// tests.cpp
#include <iostream>
#include <vector>
#include <cassert>
#include <chrono>
#include <gsl/gsl_integration.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <string>

#include "tests.hpp"
#include "profile.hpp"
#include "hydrodynamics.hpp"
#include "spectrum.hpp"
#include "maths_ops.hpp"
// #include "matplotlibcpp.h"

// namespace plt = matplotlibcpp;

/*
TO DO:
- create separate executable for tests (or maybe just run all at the start of main?)
- define size() for splines (or find other way to compare their size - probably easier)
- might not need GSL test, remove if not (and unlink in cmake)
- make more test_power robust - input x and exp and check if exp=6 to use power6
- update test_FluidProfile to test it for deflag, hybrid and detonation (i.e. repeat same tests for different PTParams input vals)
- update test_FluidProfile to take in input file to compare v,w,la profiles to xiao's code once ODE solver is working properly
*/

// Run all tests
void test_all() {
    std::cout << "Running all tests..." << "\n"
              << "Running class tests:" << std::endl;
    // test_vec();
    // test_PowerSpec();
    // test_FluidProfile();

    // finish later
    // maybe output error if individual tests fail and then continue to the next test
    // final output is either "all tests passed" or "the following tests failed: ..."

    std::cout << "All tests passed!" << std::endl;
    return;
}

// void test_rk4_coupled_odes() {
//     // dy0/dx = y1, dy1/dx = -y0 with y(0) = [1, 0] => y0(x) = cos(x), y1(x) = -sin(x)
//     auto dydx = [](double x, const state_type& y) -> state_type {
//         return { y[1], -y[0] };
//     };

//     const double x0 = 0.0;
//     const double xf = M_PI / 2;  // pi/2, where cos(pi/2) = 0, sin(pi/2) = 1
//     const size_t steps = 5000;
//     const state_type y0 = {1.0, 0.0};  // y0 = cos(0), y1 = -sin(0)

//     auto [x_vals, y_vals] = rk4_solver(dydx, x0, xf, y0, steps);

//     state_type y1_vals, y2_vals, y1_exact, y2_exact;
//     std::vector<double> err1, err2;
    
//     const auto tol = 1e-3;
//     bool test_passed = true;

//     for (int i = 0; i < x_vals.size(); i++) {
//         // fill solver result
//         y1_vals.push_back(y_vals[i][0]);
//         y2_vals.push_back(y_vals[i][1]);

//         // fill exact result
//         const auto x = x_vals[i];
//         y1_exact.push_back(std::cos(x));
//         y2_exact.push_back(-std::sin(x));

//         // error estimate
//         err1.push_back(std::abs((y1_vals[i] - y1_exact[i]) / y1_exact[i]));
//         err2.push_back(std::abs((y2_vals[i] - y2_exact[i]) / y2_exact[i]));

//         if (err1[i] > tol || err2[i] > tol) {
//             test_passed = false;
//         }
//     }

//     if (test_passed) {
//         std::cout << "Convergence test passed! (err < " << tol << ")\n";
//     } else {
//         std::cout << "Convergence test failed! (err > " << tol << ")\n";
//     }

//     // std::cout << "Plotting solution and exact result... ";

//     // plt::figure_size(2400, 600);

//     // plt::subplot2grid(1, 2, 0, 0);
//     // plt::plot(x_vals, y1_vals, "k-");
//     // plt::plot(x_vals, y1_exact, "r--");
//     // plt::xlabel("x");
//     // plt::ylabel("y0(x)");
//     // plt::grid(true);

//     // plt::subplot2grid(1, 2, 0, 1);
//     // plt::plot(x_vals, y2_vals, "k-");
//     // plt::plot(x_vals, y2_exact, "r--");
//     // plt::xlabel("x");
//     // plt::ylabel("y1(x)");
//     // plt::grid(true);

//     // const std::string filename = "test_solver.png";
//     // plt::save("../" + filename);

//     // std::cout << "Saved to file '" << filename << "'\n";
// }

void test_rk4_solver() {
    // dy/dx = -y with y(0) = 1 => y(x) = exp(-x)
    auto dydx = [](double x, const state_type& y) -> state_type {
        return {-y[0]};
    };

    const double x0 = 0.0;
    const double xf = 5.0;
    const size_t steps = 1000;
    const state_type y0 = {1.0};

    auto [x_vals, y_vals] = rk4_solver(dydx, x0, xf, y0, steps);

    // Compare final result to exact solution: y(5) = exp(-5)
    double y_numeric = y_vals.back()[0];
    double y_exact = std::exp(-xf);
    double rel_error = std::abs((y_numeric - y_exact) / y_exact);

    std::cout << "RK4 y(" << xf << ") = " << y_numeric
              << " | exact = " << y_exact
              << " | relative error = " << rel_error << '\n';

    // Assert error is small
    assert(rel_error < 1e-3);
}

// // Class tests
// void test_PowerSpec() {
//     using std::cout;
//     using std::endl;
//     using namespace Spectrum;

//     std::vector<double> k_vals{0.1, 0.2, 0.3};
//     std::vector<double> P_vals{1.0, 2.0, 3.0};

//     PowerSpec v1(k_vals, P_vals);
//     assert(v1.k().size() == 3);
//     assert(v1.P()[2] == 3.0);
//     assert(v1.max() == 3.0);

//     PowerSpec v2 = v1 * 2.0;
//     assert(v2.P()[0] == 2.0);
//     assert(v2.P()[1] == 4.0);
//     assert(v2.P()[2] == 6.0);

//     PowerSpec v3 = v1 + v1;
//     assert(v3.P()[1] == 4.0);

//     v3 *= 0.5;
//     assert(v3.P()[1] == 2.0);

//     // Error case: mismatched k-vector sizes
//     try {
//         std::vector<double> bad_k_vals{0.1, 0.2}; // shorter
//         std::vector<double> bad_P_vals{1.0, 2.0};
//         PowerSpec vbad(bad_k_vals, P_vals); // size mismatch!
//         assert(false); // Should not reach here
//     } catch (const std::invalid_argument &e) {
//         cout << "Caught expected size mismatch error: " << e.what() << endl;
//     }

//     cout << "All tests passed successfully!\n";
//     return;
// }

// void test_FluidProfile() { // not finished
//     // check read profile (from xiao's code) and solve ODE solutions are consistent
//     // plot profiles
//     using namespace Hydrodynamics;
//     using namespace PhaseTransition;

//     std::cout << "Running FluidProfile class test..." << std::endl;

//     PTParams params;
//     FluidProfile profile(params);

//     auto xi_vals = profile.xi_vals();
//     auto v_vals = profile.v_vals();
//     auto w_vals = profile.w_vals();
//     auto la_vals = profile.la_vals();

//     for (int i = 0; i < xi_vals.size(); i++) {
//         const auto xi = xi_vals[i];
//         const auto v = v_vals[i];
//         const auto w = w_vals[i];
//         const auto la = la_vals[i];

//         if (isnan(xi) || isnan(v) || isnan(w) || isnan(la)) throw std::runtime_error("NaN detected in fluid profile");

//         if (xi < 0.0) throw std::runtime_error("Domain error: xi < 0");
//         if (xi > 1.0) throw std::runtime_error("Domain error: xi > 1");

//         // also test if v,w are non-zero behind bubble wall or in front of shock
//     }

//     // profile.plot("fluid_profile_test.png");

//     std::cout << "FluidProfile test passed.\n";
// }

// update to test other interpolators? or separate one for each used?
// update to test weirder functions (or input a function to test interpolator?)
void test_interpolator() {
    // Known function: f(x) = x^3
    auto f = [](double x) { return x * x * x; };

    // Sample points
    std::vector<double> x_vals = linspace(-2.0, 2.0, 50);
    std::vector<double> y_vals;
    for (double x : x_vals)
        y_vals.push_back(f(x));

    // Build spline
    CubicSpline<double> spline(x_vals, y_vals);

    // Interpolate densely for plot
    std::vector<double> x_dense, y_interp;
    for (double x = -2.0; x <= 2.0; x += 0.01) {
        x_dense.push_back(x);
        y_interp.push_back(spline(x));
    }

    // Test interpolation accuracy at a few points
    std::vector<double> test_points = {-1.5, -0.5, 0.5, 1.5};
    double tol = 1e-4;
    for (double xi : test_points) {
        double yi_exact = f(xi);
        double yi_interp = spline(xi);
        double error = std::abs(yi_interp - yi_exact);
        std::cout << "x = " << xi
                  << ", exact = " << yi_exact
                  << ", spline = " << yi_interp
                  << ", error = " << error << '\n';
        assert(error < tol && "CubicSpline interpolation error too large!");
    }

    // Plot and save figure
    // plt::figure_size(800, 600);
    // plt::plot(x_dense, y_interp, {{"label", "Spline"}});
    // plt::scatter(x_vals, y_vals, 10.0, {{"label", "Data points"}});
    // plt::xlabel("x");
    // plt::ylabel("y");
    // plt::title("CubicSpline Interpolation of f(x) = x^3");
    // plt::legend();
    // plt::grid(true);
    // plt::save("interpolator_test.png");

    // std::cout << "CubicSpline test passed and plot saved to interpolator_test.png\n";
}

void test_prof_ints() { // test profile integrals f' and l
    PhaseTransition::PTParams params;
    Hydrodynamics::FluidProfile profile(params);

    const auto chi_vals = linspace(0.0, 100.0, 200);
    const auto [fd_vals, l_vals] = prof_ints_fl(chi_vals, profile);

    // if (plot) {
    //     plt::figure_size(1600, 600);

    //     // f'(chi)
    //     plt::subplot2grid(1, 2, 0, 0);
    //     plt::plot(chi_vals, fd_vals);
    //     plt::xlabel("chi");
    //     plt::ylabel("f'(chi)");
    //     plt::grid(true);

    //     // l(chi)
    //     plt::subplot2grid(1, 2, 0, 1);
    //     plt::plot(chi_vals, l_vals);
    //     plt::xlabel("chi");
    //     plt::ylabel("l(chi)");
    //     plt::grid(true);

    //     plt::save("fd_l_profile.png");
    // }
    return;
}

void test_Apsq() { // test Ap_sq
    PhaseTransition::PTParams params;
    Hydrodynamics::FluidProfile profile(params);

    const auto chi_vals = logspace(1e-3, 100, 1000);
    const auto Apsq = Hydrodynamics::Ap_sq(chi_vals, profile);

    // if (plot) {
    //     plt::figure_size(800, 600);
    //     plt::plot(chi_vals, Apsq);
    //     plt::xlabel("chi");
    //     plt::ylabel("Ap_sq(chi)");
    //     plt::grid(true);
    //     plt::save("Ap_sq_profile.png");
    // }
    return;
}

void test_Ekin() {
    const PhaseTransition::Universe un;

    const PhaseTransition::PTParams params1(0.5, 0.1, 1.0, 10.0, 1.71, "bag", "exp", un);
    const Hydrodynamics::FluidProfile profile1(params1);

    const PhaseTransition::PTParams params2(0.5, 0.1, 1.0, 10.0, 1.71, "bag", "sim", un);
    const Hydrodynamics::FluidProfile profile2(params2);

    const auto kRs_vals = logspace(1e-1, 1e+3, 500);

    const auto Ek1 = Spectrum::Ekin(kRs_vals, profile1);
    const auto Eks1 = Spectrum::zetaKin(Ek1);

    const auto Ek2 = Spectrum::Ekin(kRs_vals, profile2);
    const auto Eks2 = Spectrum::zetaKin(Ek2);
    
    // if (plot) {
    //     plt::figure_size(800, 600);
    //     plt::loglog(Eks1.k(), Eks1.P(), "k-");
    //     plt::loglog(Eks2.k(), Eks2.P(), "r-");
    //     plt::xlabel("K=kRs");
    //     plt::ylabel("Ekin(K)");
    //     plt::xlim(1e-1, 1e+3);
    //     plt::ylim(1e-4, 1e+0);
    //     plt::grid(true);
    //     plt::save("Ekin_spectrum.png");
    // }

    return;
}

void test_GWSpec() {
    const PhaseTransition::Universe un;

    const PhaseTransition::PTParams params1(0.5, 0.1, 1.0, 10.0, 1.71, "bag", "exp", un);
    const Hydrodynamics::FluidProfile profile1(params1);

    const PhaseTransition::PTParams params2(0.5, 0.1, 1.0, 10.0, 1.71, "bag", "sim", un);
    const Hydrodynamics::FluidProfile profile2(params2);

    const auto kRs_vec = logspace(1e-3, 1e+3, 50);
    const auto OmegaGW1 = Spectrum::GWSpec(kRs_vec, params1);
    const auto OmegaGW2 = Spectrum::GWSpec(kRs_vec, params2);
    
    // if (plot) {
    //     plt::figure_size(800, 600);
    //     plt::loglog(OmegaGW1.k(), OmegaGW1.P(), "k-");
    //     plt::loglog(OmegaGW2.k(), OmegaGW2.P(), "r-");
    //     plt::xlabel("K=kRs");
    //     plt::ylabel("Omega_GW(K)");
    //     plt::xlim(1e-3, 1e+3);
    //     // plt::ylim(1e-4, 1e+0);
    //     plt::grid(true);
    //     plt::save("GW_spectrum.png");
    // }

    return;
}

void test_SiCi() {
    // Reference values (e.g., from WolframAlpha or SciPy)
    struct TestCase {
        double x;
        double expected_si;
        double expected_ci;
    };

    TestCase tests[] = {
        {-5.0, -1.5499312449, -0.19003},
        {-4.0, -1.7582031389, -0.140982},
        {-3.0, -1.8486525274, 0.11963},
        {-2.0, -1.6054129768, 0.422981},
        {-1.0, -0.9460830704, 0.337404 },
        {-0.5, -0.4931074180, -0.177784},
        { 0.5,  0.4931074180, -0.177784},
        { 1.0,  0.9460830704, 0.337404},
        { 2.0,  1.6054129768,  0.422981},
        { 3.0,  1.8486525274,  0.11963},
        { 4.0,  1.7582031389,  -0.140982},
        { 5.0,  1.5499312449,  -0.19003}
    };

    double tol = 1e-5;

    for (const auto& t : tests) {
        // double si = Si(t.x);
        // double ci = Ci(t.x);
        const auto [si, ci] = SiCi(t.x);

        std::cout << "x = " << t.x << ": "
                  << "Si = " << si << " (expected " << t.expected_si << "), "
                  << "Ci = " << ci << " (expected " << t.expected_ci << ")\n";

        assert(std::abs(si - t.expected_si) < tol && "Si(x) does not match expected value");
        assert(std::abs(ci - t.expected_ci) < tol && "Ci(x) does not match expected value");
    }

    std::cout << "All sine/cosine integral tests passed!\n";
}

// void test_dlt_SSM() {
//     const PhaseTransition::PTParams params;

//     const auto k_vals = logspace(1e-3, 1e+3, 5);
//     const auto p_vals = linspace(1e-2, 1e+3, 200);
//     const auto z_vals = linspace(-1.0, 1.0, 200);

//     const auto nk = k_vals.size();
//     const auto np = p_vals.size();
//     const auto nz = z_vals.size();

//     const auto dlta1 = Spectrum::dlt(50, k_vals, p_vals, z_vals, params);
//     // const auto dlta2 = Spectrum::dlt_SSM(50, k_vals, p_vals, z_vals, params);

//     const auto tol = 1e-10;
//     // for (size_t kk = 0; kk < nk; kk++) {
//     //     for (size_t pp = 0; pp < np; pp++) {
//     //         for (size_t zz = 0; zz < nz; zz++) {
//     //             const auto dlt1 = dlta1[kk][pp][zz];
//     //             // const auto dlt2 = dlta2[kk][pp][zz];

//     //             // std::cout << "dlt1=" << dlt1 << ", dlt2=" << dlt2 << "\n";

//     //             // const auto error = std::abs(dlt1 - dlt2);
//     //             // if (error > tol) {
//     //             //     std::cout << "err=" << error << " for (k,p,z)=(" << kk << "," << pp << "," << zz << ")\n";
//     //             // }
//     //         }
//     //     }
//     // }

//     return;
// }

// probably don't need this
void test_simpson_integrate() {
    auto f = [](double x) { return std::sin(x); };
    auto F = [](double a, double b) { return -std::cos(b) + std::cos(a); }; // ∫ sin(x) dx

    std::cout << std::fixed << std::setprecision(10);

    // Test 1: Uniform spacing
    {
        const int N = 101;
        const double a = 0.0, b = M_PI;
        // double h = (b - a) / (N - 1);

        const auto x1_vals = linspace(a, b, N);
        std::vector<double> y1_vals;
        for (const auto x : x1_vals) {
            y1_vals.push_back(f(x));
        }

        double approx = simpson_integrate(x1_vals, y1_vals);
        double exact = F(a, b);
        std::cout << "Uniform spacing test:\n";
        std::cout << "  Simpson result: " << approx << "\n";
        std::cout << "  Exact result:   " << exact << "\n";
        std::cout << "  Error:          " << std::abs(approx - exact) << "\n";
        assert(std::abs(approx - exact) < 1e-6);
    }

    // Test 2: Non-uniform spacing
    {
        const int N = 101;
        const double a = 0.00001, b = M_PI;

        const auto x2_vals = logspace(a, b, N);
        std::vector<double> y2_vals;
        for (const auto x : x2_vals) {
            y2_vals.push_back(f(x));
        }

        double approx = simpson_integrate(x2_vals, y2_vals);
        double exact = F(a, b);
        std::cout << "\nNon-uniform spacing test:\n";
        std::cout << "  Simpson result: " << approx << "\n";
        std::cout << "  Exact result:   " << exact << "\n";
        std::cout << "  Error:          " << std::abs(approx - exact) << "\n";
        assert(std::abs(approx - exact) < 1e-4);
    }

    // Test 3: Discontinuous function (step)
    {
        auto g = [](double x) { return x < 0.5 ? 1.0 : 2.0; };
        const int N = 101;
        std::vector<double> x(N), y(N);
        for (int i = 0; i < N; ++i) {
            x[i] = static_cast<double>(i) / (N - 1);
            y[i] = g(x[i]);
        }
        double approx = simpson_integrate(x, y);
        double exact = 1.0 * 0.5 + 2.0 * 0.5; // step from 1 to 2 at x = 0.5
        std::cout << "\nDiscontinuous function test:\n";
        std::cout << "  Result: " << approx << ", Exact: " << exact << ", Error: " << std::abs(approx - exact) << "\n";
        assert(std::abs(approx - exact) < 1e-2);
    }

    // Test 5: Mismatched sizes (should throw)
    {
        std::vector<double> x = {0.0, 1.0, 2.0};
        std::vector<double> y = {1.0, 2.0}; // wrong size
        bool caught = false;
        try {
            simpson_integrate(x, y);
        } catch (const std::invalid_argument&) {
            caught = true;
        }
        assert(caught);
        std::cout << "\nMismatched size test passed (threw as expected).\n";
    }

    std::cout << "\n✅ All extended Simpson integration tests passed.\n";
}

// tests uniform/nonuniform 2d simpson integrators
void test_integrators() {
    // Known function: f(x, y) = x * y over [0,1] x [0,1]
    // Integral = ∫₀¹∫₀¹ x*y dxdy = 1/4 = 0.25

    const size_t nx = 500;
    const size_t ny = 500;
    std::vector<double> x(nx), y(ny), f_flat(nx * ny);

    // Create non-uniform x and y grids
    for (size_t i = 0; i < nx; ++i) {
        x[i] = std::pow(double(i) / (nx - 1), 1.5); // skewed toward 0
    }
    for (size_t j = 0; j < ny; ++j) {
        y[j] = std::pow(double(j) / (ny - 1), 1.2); // skewed toward 0
    }

    // Evaluate function f(x, y) = x * y on grid
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            f_flat[j * nx + i] = x[i] * y[j];
        }
    }

    double true_val = 0.25;
    double approx_uniform = simpson_2d_integrate_flat(x, y, f_flat);
    double approx_nonuniform = simpson_2d_nonuniform_flat(x, y, f_flat);

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "True integral      = " << true_val << "\n";
    std::cout << "Original integrator = " << approx_uniform << " (rel. error = " << std::abs(approx_uniform - true_val) / true_val << ")\n";
    std::cout << "Non-uniform integrator = " << approx_nonuniform << " (rel. error = " << std::abs(approx_nonuniform - true_val) / true_val << ")\n";
}
// tests for different versions of the power function
namespace PowerTest { 

double mean(const std::vector<double>& times) {
    double sum = 0.0;
    for (double time : times) {
        sum += time;
    }
    return sum / times.size();
}

double st_dev(const std::vector<double>& times, double mean_time) {
    double sum = 0.0;
    for (double time : times) {
        sum += (time - mean_time) * (time - mean_time);
    }
    return std::sqrt(sum / times.size());
}

void test_power(int num_runs) {
    double x = 2.0;  // Example input value
    int exp = 6;     // Example exponent for power()

    // Vectors to store times for each run
    std::vector<double> std_pow_times, power_times, power6_times;

    // Run std::pow multiple times and record the time
    for (int i = 0; i < num_runs; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        std::pow(x, exp);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std_pow_times.push_back(duration.count());
    }

    // Run custom power function multiple times and record the time
    for (int i = 0; i < num_runs; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        power(x, exp);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        power_times.push_back(duration.count());
    }

    // Run custom power6 function multiple times and record the time
    for (int i = 0; i < num_runs; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        power6(x);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        power6_times.push_back(duration.count());
    }

    // Calculate mean and standard deviation for each
    double mean_std_pow = mean(std_pow_times);
    double std_pow_sd = st_dev(std_pow_times, mean_std_pow);

    double mean_power = mean(power_times);
    double power_sd = st_dev(power_times, mean_power);

    double mean_power6 = mean(power6_times);
    double power6_sd = st_dev(power6_times, mean_power6);

    // Output results
    std::cout << "std::pow - Mean time: " << mean_std_pow << " seconds, SD: " << std_pow_sd << " seconds\n";
    std::cout << "Custom power function - Mean time: " << mean_power << " seconds, SD: " << power_sd << " seconds\n";
    std::cout << "Custom power6 function - Mean time: " << mean_power6 << " seconds, SD: " << power6_sd << " seconds\n";
}

} // namespace PowerTest