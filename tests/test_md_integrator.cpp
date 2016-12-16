#include <vector>

#include "tests/utils.h"

#include "src/conf.h"
#include "src/integrator.cpp"

using namespace std;

ADD_TEST(test_calc_force,
    SystemConf conf {2};
    conf.add_atom(0.0, 0.0, 0.0);
    conf.add_atom(1.0, 0.0, 0.0);

    calc_forces(conf);

    const vector<double> fs {
        1.0, 0.0, 0.0,
        -1.0, 0.0, 0.0
    };

    ASSERT_EQ_VEC(conf.fs, fs, "forces are not calculated correctly");
)

RUN_TESTS(
    test_calc_force();
)
