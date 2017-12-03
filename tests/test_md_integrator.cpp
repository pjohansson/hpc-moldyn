#include <vector>

#include "tests/utils.h"

#include "src/conf.h"
#include "src/forcefield.h"
#include "src/integrator.cpp"

using namespace std;

constexpr ForceField TestFF (
    1.0, // epsilon
    1.0, // sigma
    1.1  // rcut
);

ADD_TEST(test_calc_force,
    SystemConf conf {2};
    conf.add_atom(0.0, 0.0, 0.0);
    conf.add_atom(1.0, 0.0, 0.0);
    calc_forces_internal(conf, TestFF);


    const auto force = 24.0; // epsilon = 24, sigma = 1, dr = 1 (along x only)

    const vector<double> expected {
        force, 0.0, 0.0,
        -force, 0.0, 0.0
    };

    ASSERT_EQ_VEC(conf.fs, expected, "forces are not calculated correctly");
)

ADD_TEST(test_calc_force_outside_of_rcut_is_zero,
    SystemConf conf {2};
    conf.add_atom(0.0, 0.0, 0.0);
    conf.add_atom(TestFF.rcut + 0.1, 0.0, 0.0);

    calc_forces_internal(conf, TestFF);

    const vector<double> expected {
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    };

    ASSERT_EQ_VEC(conf.fs, expected, "force calculation does not cut at rcut");
)

ADD_TEST(test_calc_force_adds_total_force,
    SystemConf conf {2};
    conf.add_atom(0.0, 0.0, 0.0); // feels the next atom only
    conf.add_atom(1.0, 0.0, 0.0); // is pulled equally in both directions
    conf.add_atom(2.0, 0.0, 0.0); // feels the previous atom only

    calc_forces_internal(conf, TestFF);

    const auto force = 24.0;

    const vector<double> expected {
        force, 0.0, 0.0,
        0.0, 0.0, 0.0,
        -force, 0.0, 0.0
    };

    ASSERT_EQ_VEC(conf.fs, expected, "force calculation does not add correctly");
)

RUN_TESTS(
    test_calc_force();
    test_calc_force_outside_of_rcut_is_zero();
    test_calc_force_adds_total_force();
)
