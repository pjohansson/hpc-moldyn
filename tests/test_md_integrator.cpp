#include <vector>

#include "tests/utils.h"

#include "src/conf.h"
#include "src/forcefield.h"
#include "src/integrator.h"

using namespace std;

constexpr ForceField TestFF (
    1.0, // epsilon
    1.0, // sigma
    1.1, // rcut
    2.0  // mass
);

ADD_TEST(test_calc_force,
    Box box (2, RVec {0.0, 0.0, 0.0}, RVec{1.0, 1.0, 1.0});

    box.add_atom(0.0, 0.0, 0.0);
    box.add_atom(1.0, 0.0, 0.0);

    calc_forces_internal(box, TestFF);

    const auto force = 24.0; // epsilon = 24, sigma = 1, dr = 1 (along x only)

    const vector<double> expected {
        force, 0.0, 0.0,
        -force, 0.0, 0.0
    };

    ASSERT_EQ_VEC(box.fs, expected, "forces are not calculated correctly");
)

ADD_TEST(test_calc_force_outside_of_rcut_is_zero,
    Box box (2, RVec {0.0, 0.0, 0.0}, RVec{1.0, 1.0, 1.0});

    box.add_atom(0.0, 0.0, 0.0);
    box.add_atom(TestFF.rcut + 0.1, 0.0, 0.0);

    calc_forces_internal(box, TestFF);

    const vector<double> expected {
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    };

    ASSERT_EQ_VEC(box.fs, expected, "force calculation does not cut at rcut");
)

ADD_TEST(test_calc_force_adds_total_force,
    Box box (3, RVec {0.0, 0.0, 0.0}, RVec{1.0, 1.0, 1.0});

    box.add_atom(0.0, 0.0, 0.0); // feels the next atom only
    box.add_atom(1.0, 0.0, 0.0); // is pulled equally in both directions
    box.add_atom(2.0, 0.0, 0.0); // feels the previous atom only

    calc_forces_internal(box, TestFF);

    const auto force = 24.0;

    const vector<double> expected {
        force, 0.0, 0.0,
        0.0, 0.0, 0.0,
        -force, 0.0, 0.0
    };

    ASSERT_EQ_VEC(box.fs, expected, "force calculation does not add correctly");
)

ADD_TEST(test_calc_force_between_two_boxes,
    const RVec size {0.0, 0.0, 0.0};

    Box box1 (1, RVec {0.0, 0.0, 0.0}, size);
    box1.add_atom(0.0, 0.0, 0.0);

    Box box2 (1, RVec {1.0, 0.0, 0.0}, size); // shifted 1 along x
    box2.add_atom(0.0, 0.0, 0.0);

    calc_forces_from_to_box(box1, box2, TestFF);

    const auto force = 24.0; // epsilon = 24, sigma = 1, dr = 1 (along x only)

    const vector<double> expected1 {force, 0.0, 0.0};
    const vector<double> expected2 {-force, 0.0, 0.0};

    ASSERT_EQ_VEC(box1.fs, expected1, "forces are not calculated correctly");
    ASSERT_EQ_VEC(box2.fs, expected2, "forces are not calculated correctly");
)

RUN_TESTS(
    test_calc_force();
    test_calc_force_outside_of_rcut_is_zero();
    test_calc_force_adds_total_force();
    test_calc_force_between_two_boxes();
)
