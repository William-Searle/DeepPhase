#pragma once

#include <vector>
#include <matplotlibcpp.h>

void test_PowerSpec();

void test_vec();

void test_profile_solver();

double test_func(double x, void* params);

void test_boostmath();

void test_gsl_integration();

namespace PowerTest {

double mean(const std::vector<double>& times);
double st_dev(const std::vector<double>& times, double mean_time);
void test_power(int num_runs = 1000);

} // namespace PowerTest

// only takes in function of 1 variable at this stage
template<typename Func>
void plot_func(Func f, double xmin, double xmax, int points=100) {
    namespace plt = matplotlibcpp;
    
    std::vector<double> xs, ys;
    for (int i = 0; i <= points; ++i) {
        double x = xmin + i * (xmax - xmin) / points;
        xs.push_back(x);
        ys.push_back(f(x));
    }
    
    plt::plot(xs, ys);
    plt::xlabel("x");
    plt::ylabel("f(x)");
    plt::grid(true);
    // plt::show();
    plt::save("plot.png");
}

// takes input file and plots it
// void plot_pts() {

// }