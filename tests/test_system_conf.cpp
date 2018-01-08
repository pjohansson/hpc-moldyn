#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "tests/utils.h"

#include "src/conf.h"

using namespace std;

ADD_TEST(test_cell_list_init,
    const RVec origin { 0.0, 1.0, 2.0 };
    const RVec size { 3.0, 4.0, 5.0 };

    CellList list {3, origin, size};

    ASSERT_EQ(list.num_atoms(), 0, "number of atoms is not zero-initialized");
    ASSERT_EQ(list.xs.size(), 0, "number of position elements is not 0");
    ASSERT_EQ(list.vs.size(), 0, "number of velocity elements is not 0");
    ASSERT_EQ(list.fs.size(), 0, "number of force elements is not 0");
    ASSERT_EQ(list.fs_prev.size(), 0, "number of force (prev) elements is not 0");
    ASSERT_EQ(list.xs.capacity(), 9, "capacity of vector is not 3 * 3 elems as input");

    ASSERT_EQ_VEC(list.origin, origin, "list origin is not set correctly");
    ASSERT_EQ_VEC(list.size, size, "list size is not set correctly");

    ASSERT_EQ(list.to_neighbours.size(), 0, "lists are not initialized with zero neighbours as they should be");
)

ADD_TEST(test_cell_list_add_atom,
    CellList list {3, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};
    list.add_atom(1.0, 2.0, 3.0);
    list.add_atom(1.0, 2.0, 3.0);

    ASSERT_EQ(list.num_atoms(), 2, "atom number was not updated after add");
    ASSERT_EQ(list.xs.size(), 6, "wrong number of position elements");
    ASSERT_EQ(list.vs.size(), 6, "wrong number of velocity elements");
    ASSERT_EQ(list.fs.size(), 6, "wrong number of force elements");
    ASSERT_EQ(list.fs_prev.size(), 6, "wrong number of force (prev) elements");

    const vector<double> xs {1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
    const vector<double> zs {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    ASSERT_EQ_VEC(list.xs, xs, "added atom positions are not those input");
    ASSERT_EQ_VEC(list.vs, zs, "added atom velocities are not zero-initialized");
    ASSERT_EQ_VEC(list.fs, zs, "added atom forces are not zero-initialized");
    ASSERT_EQ_VEC(list.fs_prev, zs, "added atom forces (prev) are not zero-initialized");
)

ADD_TEST(test_cell_list_get_atom,
    CellList list {2, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};

    std::vector<real> xs {0.0, 1.0, 2.0, 3.0,  4.0,  5.0};
    std::vector<real> vs {0.0, 2.0, 4.0, 6.0,  8.0, 10.0};
    std::vector<real> fs {0.0, 3.0, 6.0, 9.0, 12.0, 15.0};

    list.xs = xs;
    list.vs = vs;
    list.fs = fs;

    Atom atom1 {
        RVec {0.0, 1.0, 2.0},
        RVec {0.0, 2.0, 4.0},
        RVec {0.0, 3.0, 6.0}
    };
    Atom atom2 {
        RVec {3.0,  4.0,  5.0},
        RVec {6.0,  8.0, 10.0},
        RVec {9.0, 12.0, 15.0}
    };

    ASSERT_EQ_VEC(list.get_atom(0).xs, atom1.xs, "atoms are not returned correctly");
    ASSERT_EQ_VEC(list.get_atom(0).vs, atom1.vs, "atoms are not returned correctly");
    ASSERT_EQ_VEC(list.get_atom(0).fs, atom1.fs, "atoms are not returned correctly");
    ASSERT_EQ_VEC(list.get_atom(1).xs, atom2.xs, "atoms are not returned correctly");
    ASSERT_EQ_VEC(list.get_atom(1).vs, atom2.vs, "atoms are not returned correctly");
    ASSERT_EQ_VEC(list.get_atom(1).fs, atom2.fs, "atoms are not returned correctly");
)

ADD_TEST(test_system_init,
    const std::string title {"title of system"};
    const RVec box_size {1.0, 2.0, 3.0};

    const System system {title, box_size};

    ASSERT_EQ(system.title, title, "title is not initialized correctly");
    ASSERT_EQ_VEC(system.box_size, box_size, "box size is not initialized correctly");

    ASSERT_EQ(system.num_atoms(), 0, "number of atoms is not zero-initialized");
    ASSERT_EQ(system.cell_lists.size(), 0, "number of lists in system is not zero");

    constexpr IVec shape {1, 1, 1};
    ASSERT_EQ_VEC(system.shape, shape, "shape of system is not initialized as (1, 1, 1)");
)

ADD_TEST(test_system_adds_num_atoms_from_lists,
    const std::string title {"title of system"};
    const RVec box_size {1.0, 2.0, 3.0};

    CellList list1 {1, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};
    CellList list2 {2, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};

    list1.add_atom(0.0, 0.0, 0.0);
    list1.add_atom(1.0, 0.0, 0.0);
    list2.add_atom(0.0, 0.0, 0.0);
    list2.add_atom(1.0, 0.0, 0.0);
    list2.add_atom(2.0, 0.0, 0.0);

    System system {title, box_size};
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);

    ASSERT_EQ(system.num_atoms(), 5, "system does not calculate number of atoms correctly");
)

ADD_TEST(test_system_read_grofile,
    const string path = TEST_FILES_DIRECTORY + string{"/grofile_small.gro"};
    const auto system = read_conf_from_grofile(path);

    ASSERT_EQ(system.title, "small test grofile", "incorrect title read");

    // The file has three atoms and only the positions should be read
    ASSERT_EQ(system.cell_lists.size(), 1, "the atoms were not added to a single list");
    ASSERT_EQ(system.cell_lists[0].num_atoms(), 3, "incorrect number of atoms read");

    const vector<double> xs {
        0.129, 0.079, 0.464,
        0.119, 0.061, 0.314,
        0.109, 0.042, 0.165
    };
    ASSERT_EQ_VEC(system.cell_lists[0].xs, xs, "positions were not read correctly");

    const vector<double> zeroes (9, 0.0);
    ASSERT_EQ_VEC(system.cell_lists[0].vs, zeroes, "velocities are not zero-initialized");
    ASSERT_EQ_VEC(system.cell_lists[0].fs, zeroes, "forces are not zero-initialized");
    ASSERT_EQ_VEC(system.cell_lists[0].fs_prev, zeroes, "forces (prev) are not zero-initialized");

    const array<double, 3> box_size { 216.00000, 4.67650, 110.00000 };
    ASSERT_EQ_VEC(system.box_size, box_size, "incorrect box size set to system");
    ASSERT_EQ_VEC(system.cell_lists[0].size, box_size, "incorrect box size set to cell list");
)

ADD_TEST(test_system_write_grofile,
    const string path = TEST_FILES_DIRECTORY + string{"/.out_test01.gro"};

    const std::string title {"title of system"};
    const RVec box_size {1.0, 2.0, 3.0};

    CellList list1 {1, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};
    CellList list2 {2, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};

    list1.add_atom(0.0, 1.0, 2.0);
    list1.add_atom(3.0, 4.0, 5.0);
    list2.add_atom(6.0, 7.0, 8.0);

    System output {title, box_size};
    output.cell_lists.push_back(list1);
    output.cell_lists.push_back(list2);

    write_conf_to_grofile(output, path);
    const auto system = read_conf_from_grofile(path);

    ASSERT_EQ(system.title, title, "title was not written correctly");
    ASSERT_EQ_VEC(system.box_size, box_size, "box dimensions were not written correctly");

    const vector<real> xs {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    ASSERT_EQ_VEC(system.cell_lists[0].xs, xs, "positions were not written correctly");
)

ADD_TEST(test_split_system_into_lists_of_size_2x_rcut_gives_correct_shape,
    const std::string title {"title of system"};
    const RVec box_size {10.0, 20.0, 30.0};
    const float rcut = 1.1;

    const auto nx = static_cast<uint64_t>(floor(box_size[XX] / (2.0 * rcut)));
    const auto ny = static_cast<uint64_t>(floor(box_size[YY] / (2.0 * rcut)));
    const auto nz = static_cast<uint64_t>(floor(box_size[ZZ] / (2.0 * rcut)));
    const auto expected_shape = IVec {nx, ny, nz};

    const auto dx = box_size[XX] / nx;
    const auto dy = box_size[YY] / ny;
    const auto dz = box_size[ZZ] / nz;
    const auto small_box_size = RVec {dx, dy, dz};

    auto system = System(title, box_size);
    create_cell_lists(system, rcut);

    ASSERT_EQ_VEC(system.shape, expected_shape, "system does not split into the correct shape");
    ASSERT_EQ(nx * ny * nz, system.cell_lists.size(), "system does not split into the correct number of cell lists");

    for (unsigned ix = 0; ix < nx; ++ix)
    {
        for (unsigned iy = 0; iy < ny; ++iy)
        {
            for (unsigned iz = 0; iz < nz; ++iz)
            {
                const auto i = ix * ny * nz + iy * nz + iz;
                const auto x0 = ix * dx;
                const auto y0 = iy * dy;
                const auto z0 = iz * dz;
                const auto origin = RVec {x0, y0, z0};

                ASSERT_EQ_VEC(origin, system.cell_lists[i].origin,
                    "the system cell lists are not correctly ordered or "
                    "the origin is not set correctly in all of them");
                ASSERT_EQ_VEC(small_box_size, system.cell_lists[i].size,
                    "the system cell lists are not set to the correct box size");
            }
        }
    }
)

ADD_TEST(test_split_system_into_lists_creates_minimum_one_per_size,
    const std::string title {"title of system"};
    const RVec box_size {1.0, 1.0, 1.0};
    const float rcut = 1.1;

    auto system = System(title, box_size);
    create_cell_lists(system, rcut);

    const auto shape = IVec {1, 1, 1};
    ASSERT_EQ_VEC(shape, system.shape,
        "system does not split into a minimum of one cell per side "
        "which it should");
    ASSERT_EQ(1, system.cell_lists.size(),
        "system does not split into a minimum of one cell per side "
        "which it should");
    ASSERT_EQ_VEC(box_size, system.cell_lists[0].size,
        "splitting into a single cell per side gives incorrect cell sizes");

    const auto origin = RVec {0.0, 0.0, 0.0};
    ASSERT_EQ_VEC(origin, system.cell_lists[0].origin,
        "splitting into a single cell per side gives incorrect origin");
)

ADD_TEST(test_split_system_puts_atoms_in_correct_lists,
    const std::string title {"title of system"};
    const RVec box_size {4.1, 1.0, 1.0}; // split into (2, 1, 1)
    const float rcut = 1.0;

    const auto expected_shape = IVec {2, 1, 1};

    const auto dx = box_size[XX] / 2;
    const auto dy = box_size[YY];
    const auto dz = box_size[ZZ];
    const auto small_box_size = RVec {dx, dy, dz};

    // The system has three lists (of the same size, just to ensure that
    // they are properly collapsed) with one atom each, one of which should
    // be in the second cell of the final system.

    auto list1 = CellList(1, RVec {0.0, 0.0, 0.0}, box_size);
    list1.add_atom(dx / 2.0, dy / 2.0, dy / 2.0); // Should be in the first list

    auto list2 = CellList(1, RVec {0.0, 0.0, 0.0}, box_size);
    list2.add_atom(dx / 3.0, dy / 2.0, dy / 2.0); // Should be in the first list

    auto list3 = CellList(1, RVec {0.0, 0.0, 0.0}, box_size);
    list3.add_atom(1.5 * dx, dy / 2.0, dy / 2.0); // Should be in the second list

    auto system = System(title, box_size);
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);
    system.cell_lists.push_back(list3);

    create_cell_lists(system, rcut);

    ASSERT_EQ_VEC(system.shape, expected_shape,
        "this tests hard coded parameters do not give the expected shape, "
        "fix *the test*");

    ASSERT_EQ(system.cell_lists[0].num_atoms(), 2,
        "the first cell list does not contain the correct number of atoms");
    ASSERT_EQ(system.cell_lists[1].num_atoms(), 1,
        "the second cell list does not contain the correct number of atoms");

    // First list is at (0, 0, 0) so these coordinates do not have to be shifted.
    const auto xs1 = vector<real> {
        dx / 2.0, dy / 2.0, dz / 2.0,
        dx / 3.0, dy / 2.0, dz / 2.0
    };
    const auto zeroes1 = vector<real>(6, 0.0);
    ASSERT_EQ_VEC(system.cell_lists[0].xs, xs1,
        "atom coordinates are not set correctly when splitting");
    ASSERT_EQ_VEC(system.cell_lists[0].vs, zeroes1,
        "atom velocities are not set correctly when splitting");
    ASSERT_EQ_VEC(system.cell_lists[0].fs, zeroes1,
        "atom forces are not set correctly when splitting");
    ASSERT_EQ_VEC(system.cell_lists[0].fs_prev, zeroes1,
        "atom forces (prev) are not set correctly when splitting");

    // Second list is at (small_box_size[XX], 0, 0) so we need to adjust
    // the x coordinate from that in absolute space to relative
    const auto xs2 = vector<real> {1.5 * dx - small_box_size[XX], dy / 2.0, dz / 2.0};
    const auto zeroes2 = vector<real>(3, 0.0);
    ASSERT_EQ_VEC(system.cell_lists[1].xs, xs2,
        "atom coordinates are not set correctly when splitting");
    ASSERT_EQ_VEC(system.cell_lists[1].vs, zeroes2,
        "atom velocities are not set correctly when splitting");
    ASSERT_EQ_VEC(system.cell_lists[1].fs, zeroes2,
        "atom forces are not set correctly when splitting");
    ASSERT_EQ_VEC(system.cell_lists[1].fs_prev, zeroes2,
        "atom forces (prev) are not set correctly when splitting");

    // Ensure that we use the correct amount of memory
    ASSERT_EQ(system.cell_lists[0].xs.size(), system.cell_lists[0].xs.capacity(),
        "after readding the atoms the reserved memory has not been minimized");
    ASSERT_EQ(system.cell_lists[0].vs.size(), system.cell_lists[0].vs.capacity(),
        "after readding the atoms the reserved memory has not been minimized");
    ASSERT_EQ(system.cell_lists[0].fs.size(), system.cell_lists[0].fs.capacity(),
        "after readding the atoms the reserved memory has not been minimized");
    ASSERT_EQ(system.cell_lists[0].fs_prev.size(), system.cell_lists[0].fs_prev.capacity(),
        "after readding the atoms the reserved memory has not been minimized");
    ASSERT_EQ(system.cell_lists[1].xs.size(), system.cell_lists[1].xs.capacity(),
        "after readding the atoms the reserved memory has not been minimized");
    ASSERT_EQ(system.cell_lists[1].vs.size(), system.cell_lists[1].vs.capacity(),
        "after readding the atoms the reserved memory has not been minimized");
    ASSERT_EQ(system.cell_lists[1].fs.size(), system.cell_lists[1].fs.capacity(),
        "after readding the atoms the reserved memory has not been minimized");
    ASSERT_EQ(system.cell_lists[1].fs_prev.size(), system.cell_lists[1].fs_prev.capacity(),
        "after readding the atoms the reserved memory has not been minimized");
)

RUN_TESTS(
    test_cell_list_init();
    test_cell_list_add_atom();
    test_cell_list_get_atom();
    test_system_init();
    test_system_adds_num_atoms_from_lists();
    test_system_read_grofile();
    test_system_write_grofile();
    test_split_system_into_lists_of_size_2x_rcut_gives_correct_shape();
    test_split_system_into_lists_creates_minimum_one_per_size();
    test_split_system_puts_atoms_in_correct_lists();
);
