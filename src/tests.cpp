#include <iostream>
#include <vector>
#include <cassert>

#include "spectrum.hpp"
#include "spectrum.tpp"

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
}