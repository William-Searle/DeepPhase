#include "catch/catch.hpp"
#include "spectrum.hpp"

TEST_CASE("Tests functionality of the PowerSpec class.", "[powerSpec]") {
    using namespace Spectrum;

    std::vector<double> k_vals{0.1, 0.2, 0.3};
    std::vector<double> P_vals{1.0, 2.0, 3.0};

    PowerSpec v1(k_vals, P_vals);
    REQUIRE(v1.k().size() == 3);
    REQUIRE(v1.P()[2] == 3.0);
    REQUIRE(v1.max() == 3.0);

    PowerSpec v2 = v1 * 2.0;
    REQUIRE(v2.P()[0] == 2.0);
    REQUIRE(v2.P()[1] == 4.0);
    REQUIRE(v2.P()[2] == 6.0);

    PowerSpec v3 = v1 + v1;
    REQUIRE(v3.P()[1] == 4.0);

    v3 *= 0.5;
    REQUIRE(v3.P()[1] == 2.0);

    // Error case: mismatched k-vector sizes
    REQUIRE_THROWS_AS(PowerSpec({0.1, 0.2}, P_vals), std::invalid_argument);
}