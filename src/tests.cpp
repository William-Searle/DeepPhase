#include <iostream>
#include <vector>
#include <cassert>
#include <chrono>
#include <gsl/gsl_integration.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "spectrum.hpp"
#include "maths_ops.hpp"
#include "rk4_solver.hpp"

/*
TO DO:
- create separate executable for tests (i.e. have a second main function)
*/

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

// void test_profile_solver() {
//     RK4::State init{1.0, 0.0, 0.0}; // Initial condition: x=1, v=0, w=0
//     double t0 = 0.0;
//     double tf = 1.0;
//     int n = 100;

//     const auto prof = RK4::profile(RK4::State& init, double t0, double tf, int n);
//     std::cout << "Final state:\n";
//     std::cout << "xi = " << prof.xi << "\n";
//     std::cout << "v = " << prof.v << "\n";
//     std::cout << "w = " << prof.w << "\n";
// }

double test_func(double x, void *params) {
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
    gsl_func.function = &test_func;  // Assign the function pointer
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

namespace PowerTest { // test for different versions of the power function

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

// TO DO: make more robust - input x and exp and check if exp=6 to use power6
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