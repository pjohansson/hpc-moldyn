#include <array>
#include <string>
#include <vector>

#include "tests/utils.h"

#include "src/conf.h"

using namespace std;

ADD_TEST(test_box_init,
    const RVec origin { 0.0, 1.0, 2.0 };
    const RVec size { 3.0, 4.0, 5.0 };

    Box box {3, origin, size};

    ASSERT_EQ(box.num_atoms(), 0, "number of atoms is not zero-initialized");
    ASSERT_EQ(box.xs.size(), 0, "number of position elements is not 0");
    ASSERT_EQ(box.vs.size(), 0, "number of velocity elements is not 0");
    ASSERT_EQ(box.fs.size(), 0, "number of force elements is not 0");
    ASSERT_EQ(box.fs_prev.size(), 0, "number of force (prev) elements is not 0");
    ASSERT_EQ(box.xs.capacity(), 9, "capacity of vector is not 3 * 3 elems as input");

    ASSERT_EQ_VEC(box.origin, origin, "box origin is not set correctly");
    ASSERT_EQ_VEC(box.size, size, "box size is not set correctly");

    ASSERT_EQ(box.to_neighbours.size(), 0, "boxes are not initialized with zero neighbours as they should be");
)

ADD_TEST(test_box_add_atom,
    Box box {3, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};
    box.add_atom(1.0, 2.0, 3.0);
    box.add_atom(1.0, 2.0, 3.0);

    ASSERT_EQ(box.num_atoms(), 2, "atom number was not updated after add");
    ASSERT_EQ(box.xs.size(), 6, "wrong number of position elements");
    ASSERT_EQ(box.vs.size(), 6, "wrong number of velocity elements");
    ASSERT_EQ(box.fs.size(), 6, "wrong number of force elements");
    ASSERT_EQ(box.fs_prev.size(), 6, "wrong number of force (prev) elements");

    const vector<double> xs {1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
    const vector<double> zs {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    ASSERT_EQ_VEC(box.xs, xs, "added atom positions are not those input");
    ASSERT_EQ_VEC(box.vs, zs, "added atom velocities are not zero-initialized");
    ASSERT_EQ_VEC(box.fs, zs, "added atom forces are not zero-initialized");
    ASSERT_EQ_VEC(box.fs_prev, zs, "added atom forces (prev) are not zero-initialized");
)

ADD_TEST(test_system_init,
    const std::string title {"title of system"};
    const RVec box_size {1.0, 2.0, 3.0};

    const System system {title, box_size};

    ASSERT_EQ(system.title, title, "title is not initialized correctly");
    ASSERT_EQ_VEC(system.box_size, box_size, "box size is not initialized correctly");

    ASSERT_EQ(system.num_atoms(), 0, "number of atoms is not zero-initialized");
    ASSERT_EQ(system.boxes.size(), 0, "number of boxes in system is not zero");

    constexpr IVec shape {1, 1, 1};
    ASSERT_EQ_VEC(system.shape, shape, "shape of system is not initialized as (1, 1, 1)");
)

ADD_TEST(test_system_adds_num_atoms_from_boxes,
    const std::string title {"title of system"};
    const RVec box_size {1.0, 2.0, 3.0};

    Box box1 {1, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};
    Box box2 {2, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};

    box1.add_atom(0.0, 0.0, 0.0);
    box1.add_atom(1.0, 0.0, 0.0);
    box2.add_atom(0.0, 0.0, 0.0);
    box2.add_atom(1.0, 0.0, 0.0);
    box2.add_atom(2.0, 0.0, 0.0);

    System system {title, box_size};
    system.boxes.push_back(box1);
    system.boxes.push_back(box2);

    ASSERT_EQ(system.num_atoms(), 5, "system does not calculate number of atoms correctly");
)

ADD_TEST(test_system_read_grofile,
    const string path = TEST_FILES_DIRECTORY + string{"/grofile_small.gro"};
    const auto system = read_conf_from_grofile(path);

    ASSERT_EQ(system.title, "small test grofile", "incorrect title read");

    // The file has three atoms and only the positions should be read
    ASSERT_EQ(system.boxes.size(), 1, "the atoms were not added to a single box");
    ASSERT_EQ(system.boxes[0].num_atoms(), 3, "incorrect number of atoms read");

    const vector<double> xs {
        0.129, 0.079, 0.464,
        0.119, 0.061, 0.314,
        0.109, 0.042, 0.165
    };
    ASSERT_EQ_VEC(system.boxes[0].xs, xs, "positions were not read correctly");

    const vector<double> zeroes (9, 0.0);
    ASSERT_EQ_VEC(system.boxes[0].vs, zeroes, "velocities are not zero-initialized");
    ASSERT_EQ_VEC(system.boxes[0].fs, zeroes, "forces are not zero-initialized");
    ASSERT_EQ_VEC(system.boxes[0].fs_prev, zeroes, "forces (prev) are not zero-initialized");

    const array<double, 3> box_size { 216.00000, 4.67650, 110.00000 };
    ASSERT_EQ_VEC(system.box_size, box_size, "incorrect box size set to system");
    ASSERT_EQ_VEC(system.boxes[0].size, box_size, "incorrect box size set to box");
)

ADD_TEST(test_system_write_grofile,
    const string path = TEST_FILES_DIRECTORY + string{"/.out_test01.gro"};

    const std::string title {"title of system"};
    const RVec box_size {1.0, 2.0, 3.0};

    Box box1 {1, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};
    Box box2 {2, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};

    box1.add_atom(0.0, 1.0, 2.0);
    box1.add_atom(3.0, 4.0, 5.0);
    box2.add_atom(6.0, 7.0, 8.0);

    System output {title, box_size};
    output.boxes.push_back(box1);
    output.boxes.push_back(box2);

    write_conf_to_grofile(output, path);
    const auto system = read_conf_from_grofile(path);

    ASSERT_EQ(system.title, title, "title was not written correctly");
    ASSERT_EQ_VEC(system.box_size, box_size, "box dimensions were not written correctly");

    const vector<real> xs {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    ASSERT_EQ_VEC(system.boxes[0].xs, xs, "positions were not written correctly");
)

RUN_TESTS(
    test_box_init();
    test_box_add_atom();
    test_system_init();
    test_system_adds_num_atoms_from_boxes();
    test_system_read_grofile();
    test_system_write_grofile();
);
