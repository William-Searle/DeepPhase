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
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

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
    test_vec();
    test_PowerSpec();
    test_FluidProfile();

    // finish later
    // maybe output error if individual tests fail and then continue to the next test
    // final output is either "all tests passed" or "the following tests failed: ..."

    std::cout << "All tests passed!" << std::endl;
    return;
}

// Class tests
void test_PowerSpec() {
    using std::cout;
    using std::endl;
    using namespace Spectrum;

    // Scalar tests
    PowerSpec s1(1.0, 10.0);
    PowerSpec s2(1.0, 5.0);

    assert(s1.k() == 1.0);
    assert(s1.P() == 10.0);
    assert(s1.max() == 10.0);
    assert(s1.is_scalar());

    // Scalar arithmetic
    PowerSpec s3 = s1 + 2.0;
    assert(s3.P() == 12.0);

    s3 -= 2.0;
    assert(s3.P() == 10.0);

    PowerSpec s4 = s1 + s2;
    assert(s4.P() == 15.0);

    // Vector tests
    std::vector<double> kvec{0.1, 0.2, 0.3};
    std::vector<double> Pvec{1.0, 2.0, 3.0};

    PowerSpec v1(kvec, Pvec);
    assert(!v1.is_scalar());
    assert(v1.kvec().size() == 3);
    assert(v1.Pvec()[2] == 3.0);
    assert(v1.max() == 3.0);

    PowerSpec v2 = v1 * 2.0;
    assert(v2.Pvec()[0] == 2.0);
    assert(v2.Pvec()[1] == 4.0);
    assert(v2.Pvec()[2] == 6.0);

    PowerSpec v3 = v1 + v1;
    assert(v3.Pvec()[1] == 4.0);

    v3 *= 0.5;
    assert(v3.Pvec()[1] == 2.0);

    // Error case: mismatched k-vector sizes
    try {
        std::vector<double> bad_kvec{0.1, 0.2}; // shorter
        std::vector<double> bad_Pvec{1.0, 2.0};
        PowerSpec vbad(bad_kvec, Pvec); // size mismatch!
        assert(false); // Should not reach here
    } catch (const std::invalid_argument &e) {
        cout << "Caught expected size mismatch error: " << e.what() << endl;
    }

    // Error case: mixing scalar and vector
    try {
        PowerSpec bad = s1 + v1;
        assert(false); // Should not reach here
    } catch (const std::invalid_argument &e) {
        cout << "Caught expected scalar/vector mix error: " << e.what() << endl;
    }

    cout << "All tests passed successfully!\n";
    return;
}

void test_vec() {
    vec<double> a{1.0, 2.0, 3.0};
    vec<double> b{4.0, 5.0, 6.0};

    auto c = a + b;
    auto d = a * 2.0;
    auto dot = a.dot(b);

    c.print();                           // [ 5 7 9 ]
    d.print();                           // [ 2 4 6 ]
    std::cout << "Dot: " << dot << "\n"; // Dot: 32

    return;
}

void test_FluidProfile() { // not finished
    // check read profile (from xiao's code) and solve ODE solutions are consistent
    // plot profiles
    using namespace Hydrodynamics;
    using namespace PhaseTransition;

    std::cout << "Running FluidProfile class test..." << std::endl;

    PTParams params;
    FluidProfile profile(params);

    // Generate streamplot data
    try {
        profile.generate_streamplot_data(30, 30, "test_streamplot.csv");
    } catch (...) {
        assert(false && "generate_streamplot_data threw an exception");
    }

    auto xi_vals = profile.xi_vals();
    auto v_interp = profile.v_prof();
    auto w_interp = profile.w_prof();
    auto la_interp = profile.la_prof();

    profile.plot("fluid_profile_test.png");

    std::cout << "FluidProfile test passed.\n";
}

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
    plt::figure_size(800, 600);
    plt::plot(x_dense, y_interp, {{"label", "Spline"}});
    plt::scatter(x_vals, y_vals, 10.0, {{"label", "Data points"}});
    plt::xlabel("x");
    plt::ylabel("y");
    plt::title("CubicSpline Interpolation of f(x) = x^3");
    plt::legend();
    plt::grid(true);
    plt::save("interpolator_test.png");

    std::cout << "CubicSpline test passed and plot saved to interpolator_test.png\n";
}

// Integration tests
double test_func_gsl(double x, void *params) {
    // integral over -inf to inf gives sqrt(pi) = 1.77245
    return std::exp(-x * x);
}

void test_boostmath() {
    auto f = [](double x) { return std::exp(-x*x); };

    boost::math::quadrature::gauss_kronrod<double, 15> intgrl_type;
    double result = intgrl_type.integrate(f, -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
    std::cout << "result: " << result << "\n";
    return;
}

void test_gsl_integration() { // might not need GSL, remove if not (unlink in cmake)
    // Create a GSL integration workspace
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

    // Define the result variables
    double result, error;

    // Wrap the function in a gsl_function struct
    gsl_function gsl_func;
    gsl_func.function = &test_func_gsl;  // Assign the function pointer
    gsl_func.params = nullptr;       // Pass any additional parameters (none in this case)

    // Perform the integration using the GSL integration function
    gsl_integration_qagi(&gsl_func, 0, 1e-7, 1000, workspace, &result, &error);

    // Output the result and the estimated error
    std::cout << "Result of the integration: " << result << std::endl;
    std::cout << "Estimated error: " << error << std::endl;

    // Free the GSL workspace
    gsl_integration_workspace_free(workspace);

    return;
}

void test_prof_ints(bool plot) { // test profile integrals f' and l
    PhaseTransition::PTParams params;
    Hydrodynamics::FluidProfile profile(params);

    const auto chi_vals = linspace(0.0, 100.0, 200);
    const auto [fd_vals, l_vals] = prof_ints_fl(chi_vals, profile);

    if (plot) {
        plt::figure_size(1600, 600);

        // f'(chi)
        plt::subplot2grid(1, 2, 0, 0);
        plt::plot(chi_vals, fd_vals);
        plt::xlabel("chi");
        plt::ylabel("f'(chi)");
        plt::grid(true);

        // l(chi)
        plt::subplot2grid(1, 2, 0, 1);
        plt::plot(chi_vals, l_vals);
        plt::xlabel("chi");
        plt::ylabel("l(chi)");
        plt::grid(true);

        plt::save("fd_l_profile.png");
    }
    return;
}

void test_Apsq(bool plot) { // test Ap_sq
    PhaseTransition::PTParams params;
    Hydrodynamics::FluidProfile profile(params);

    // const auto chi_vals = logspace(0.1, 10, 1000);
    const auto chi_vals = linspace(0.0, 100.0, 200);
    const auto Apsq = Hydrodynamics::Ap_sq(chi_vals, profile);

    if (plot) {
        plt::figure_size(800, 600);
        plt::plot(chi_vals, Apsq);
        plt::xlabel("chi");
        plt::ylabel("Ap_sq(chi)");
        plt::grid(true);
        plt::save("Ap_sq_profile.png");
    }
    return;
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
        double result_std_pow = std::pow(x, exp);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std_pow_times.push_back(duration.count());
    }

    // Run custom power function multiple times and record the time
    for (int i = 0; i < num_runs; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        double result_power = power(x, exp);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        power_times.push_back(duration.count());
    }

    // Run custom power6 function multiple times and record the time
    for (int i = 0; i < num_runs; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        double result_power6 = power6(x);
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