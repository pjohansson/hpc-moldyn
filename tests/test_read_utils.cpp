#include "tests/utils.h"

#include <string>

#include "src/params.cpp"

using namespace std;

ADD_TEST(test_read_full_parameter_file,
    const string path = TEST_FILES_DIRECTORY + string { "/params.dat" };

    Options opts;
    ASSERT(read_parameter_file(path, opts),
        "parameter file could not be opened");

    ASSERT_EQ(opts.dt, 0.005, "dt was not read correctly");
    ASSERT_EQ(opts.energy_calc, 40, "energy_calc was not read correctly");
    ASSERT_EQ(opts.traj_stride, 60, "traj_stride was not read correctly");
    ASSERT_EQ(opts.num_steps, 2000, "num_steps was not read correctly");
    ASSERT_EQ(opts.gen_velocities, true,
        "gen_velocities was not read correctly");
    ASSERT_EQ(opts.gen_temp, 20.0, "gen_temp was not read correctly");

    // Should not be set by the file
    ASSERT_EQ(opts.verbose, false, "verbose was not initialized to false");
)

ADD_TEST(test_read_bad_parameter_file_catches_exception,
    const string path = TEST_FILES_DIRECTORY + string { "/params_bad.dat" };
    Options opts;
    ASSERT(!read_parameter_file(path, opts), "an error was not returned");
)

ADD_TEST(test_read_non_existing_file_catches_exception,
    const string path = TEST_FILES_DIRECTORY + string { "/does_not_exist" };
    Options opts;
    ASSERT(!read_parameter_file(path, opts), "an error was not returned");
)

RUN_TESTS(
    test_read_full_parameter_file();
    test_read_bad_parameter_file_catches_exception();
    test_read_non_existing_file_catches_exception();
)
