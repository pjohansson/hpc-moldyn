#include <string>

#include "tests/utils.h"

#include "src/conf.h"

using namespace std;

void test_system_init()
{
    SystemConf conf {3};
    ASSERT_EQ(conf.num_atoms(), 0, "number of atoms is not zero-initialized");
    ASSERT_EQ(conf.xs.size(), 0, "number of position elements is not 0");
    ASSERT_EQ(conf.vs.size(), 0, "number of velocity elements is not 0");
    ASSERT_EQ(conf.fs.size(), 0, "number of force elements is not 0");
    ASSERT_EQ(conf.xs.capacity(), 9, "capacity of vector is not 3*3 elems as input");

    const array<double, 3> box { 0.0, 0.0, 0.0 };
    ASSERT_EQ_VEC(conf.box, box, "box is not zero-initialized");

    ASSERT_EQ(conf.title, "", "title is not initialized as empty");

    finalize_test();
}

void test_system_add_atom()
{
    SystemConf conf {3};
    auto num = conf.add_atom(1.0, 2.0, 3.0);
    ASSERT_EQ(num, 1, "new atom number was not correctly returned after add");
    num = conf.add_atom(1.0, 2.0, 3.0);
    ASSERT_EQ(num, 2, "new atom number was not correctly returned after add");

    ASSERT_EQ(conf.num_atoms(), 2, "atom number was not updated after add");
    ASSERT_EQ(conf.xs.size(), 6, "wrong number of position elements");
    ASSERT_EQ(conf.vs.size(), 6, "wrong number of velocity elements");
    ASSERT_EQ(conf.fs.size(), 6, "wrong number of force elements");

    const vector<double> xs {1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
    const vector<double> zs {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    ASSERT_EQ_VEC(conf.xs, xs, "added atom positions are not those input");
    ASSERT_EQ_VEC(conf.vs, zs, "added atom velocities are not zero-initialized");
    ASSERT_EQ_VEC(conf.fs, zs, "added atom forces are not zero-initialized");

    finalize_test();
}

void test_system_set_box()
{
    SystemConf conf {1};
    conf.set_box(1.0, 2.0, 3.0);
    const array<double, 3> box {1.0, 2.0, 3.0};
    ASSERT_EQ_VEC(conf.box, box, "box was not set correctly");

    finalize_test();
}

void test_system_read_grofile()
{
    const string path = TEST_FILES_DIRECTORY + string{"/grofile_small.gro"};
    const auto conf = read_conf_from_grofile(path);

    ASSERT_EQ(conf.title, "small test grofile", "incorrect title read");

    // The file has three atoms and only the positions should be read
    ASSERT_EQ(conf.num_atoms(), 3, "incorrect number of atoms read");

    const vector<double> xs {
        0.129, 0.079, 0.464,
        0.119, 0.061, 0.314,
        0.109, 0.042, 0.165
    };
    ASSERT_EQ_VEC(conf.xs, xs, "positions were not read correctly");

    const vector<double> zs (9, 0.0);
    ASSERT_EQ_VEC(conf.vs, zs, "velocities are not zero-initialized");
    ASSERT_EQ_VEC(conf.fs, zs, "forces are not zero-initialized");

    const array<double, 3> box { 216.00000, 4.67650, 110.00000 };
    ASSERT_EQ_VEC(conf.box, box, "incorrect box size read");

    finalize_test();
}

RUN_TESTS(
    test_system_init();
    test_system_add_atom();
    test_system_set_box();
    test_system_read_grofile();
);
