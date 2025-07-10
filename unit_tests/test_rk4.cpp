#include "catch/catch.hpp"
#include "maths_ops.hpp"

TEST_CASE("Test rk4 Coupled ODEs", "[rk4]") {

    // dy0/dx = y1, dy1/dx = -y0 with y(0) = [1, 0] => y0(x) = cos(x), y1(x) = -sin(x)

    auto dydx = [](double x, const state_type& y) -> state_type {
        return { y[1], -y[0] };
    };

    const double x0 = 0.0;
    const double xf = M_PI / 2; // pi/2
    const size_t steps = 5000;
    const double tol = 1e-3;
    const double epsilon = 1e-10; // for safe relative error denom

    const state_type y0 = {1.0, 0.0};

    auto [x_vals, y_vals] = rk4_solver(dydx, x0, xf, y0, steps);

    for (size_t i = 0; i < x_vals.size(); ++i) {
        double x = x_vals[i];
        double y0_num = y_vals[i][0];
        double y1_num = y_vals[i][1];

        double y0_exact = std::cos(x);
        double y1_exact = -std::sin(x);

        double denom0 = std::abs(y0_exact) > epsilon ? std::abs(y0_exact) : 1.0;
        double denom1 = std::abs(y1_exact) > epsilon ? std::abs(y1_exact) : 1.0;

        double rel_err0 = std::abs(y0_num - y0_exact) / denom0;
        double rel_err1 = std::abs(y1_num - y1_exact) / denom1;

        REQUIRE(rel_err0 <= tol);
        REQUIRE(rel_err1 <= tol);
    }

}

TEST_CASE("Test rk4 Solver", "[rk4]") {

    // dy0/dx = y1, dy1/dx = -y0 with y(0) = [1, 0] => y0(x) = cos(x), y1(x) = -sin(x)

    auto dydx = [](double x, const state_type& y) -> state_type {
        return { y[1], -y[0] };
    };

    const double x0 = 0.0;
    const double xf = M_PI / 2; // pi/2
    const size_t steps = 5000;
    const double tol = 1e-3;
    const double epsilon = 1e-10; // for safe relative error denom

    const state_type y0 = {1.0, 0.0};

    auto [x_vals, y_vals] = rk4_solver(dydx, x0, xf, y0, steps);

    for (size_t i = 0; i < x_vals.size(); ++i) {
        double x = x_vals[i];
        double y0_num = y_vals[i][0];
        double y1_num = y_vals[i][1];

        double y0_exact = std::cos(x);
        double y1_exact = -std::sin(x);

        double denom0 = std::abs(y0_exact) > epsilon ? std::abs(y0_exact) : 1.0;
        double denom1 = std::abs(y1_exact) > epsilon ? std::abs(y1_exact) : 1.0;

        double rel_err0 = std::abs(y0_num - y0_exact) / denom0;
        double rel_err1 = std::abs(y1_num - y1_exact) / denom1;

        REQUIRE(rel_err0 <= tol);
        REQUIRE(rel_err1 <= tol);
    }
    
}