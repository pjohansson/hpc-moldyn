#include <vector>

#include "tests/utils.h"

#include "src/analytics.cpp"
#include "src/conf.h"
#include "src/params.h"
#include "src/integrator.cpp"

using namespace std;

constexpr ForceField TestFF (
    1.0, // epsilon
    1.0, // sigma
    1.1, // rcut
    2.0, // mass
    5.0 // wall_constant
);

const Options TestOpts;

static MPIRank get_single_rank_mpi_comm(const System& system)
{
    const auto num_cells = system.cell_lists.size();

    MPIRank mpi_comm { 0, 1 };

    vector<vector<size_t>> owned_cells (1);

    for (size_t i = 0; i < num_cells; ++i)
    {
        owned_cells.at(0).push_back(i);
    }

    // vector<vector_size_t>> non_owned_cells (1);

    mpi_comm.mpi_rank_owned_cells = owned_cells;
    // mpi_comm.mpi_rank

    return mpi_comm;
}

ADD_TEST(test_calc_force,
    constexpr ForceField TestFF_long_rcut (
        1.0, // epsilon
        1.0, // sigma
        10.0, // rcut
        2.0, // mass
        5.0 // wall_constant
    );

    CellList list (2, RVec {0.0, 0.0, 0.0}, RVec{10.0, 1.0, 1.0});

    const auto dx = 5.0;
    list.add_atom(0.0, 0.0, 0.0);
    list.add_atom(dx, 0.0, 0.0);

    calc_forces_internal(list, TestFF_long_rcut);

    const auto force = -48.0 * dx
        * (1.0 / std::pow(dx, 14) - 0.5 / std::pow(dx, 8));

    const vector<double> expected {
        force, 0.0, 0.0,
        -force, 0.0, 0.0
    };

    ASSERT_EQ_VEC(list.fs, expected, "forces are not calculated correctly");
)

ADD_TEST(test_calc_force_with_no_atoms_in_cell_works,
    // Ensure that the internal loop calculation works for the edge
    // case of no atoms in the cell.
    CellList list (0, RVec {0.0, 0.0, 0.0}, RVec{1.0, 1.0, 1.0});

    calc_forces_internal(list, TestFF);

    const vector<double> expected;

    ASSERT_EQ_VEC(list.fs, expected, "forces are not calculated correctly");
)

ADD_TEST(test_calc_force_outside_of_rcut_is_zero,
    CellList list (2, RVec {0.0, 0.0, 0.0}, RVec{1.0, 1.0, 1.0});

    list.add_atom(0.0, 0.0, 0.0);
    list.add_atom(TestFF.rcut + 0.1, 0.0, 0.0);

    calc_forces_internal(list, TestFF);

    const vector<double> expected {
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    };

    ASSERT_EQ_VEC(list.fs, expected, "force calculation does not cut at rcut");
)

ADD_TEST(test_calc_force_adds_total_force,
    CellList list (3, RVec {0.0, 0.0, 0.0}, RVec{1.0, 1.0, 1.0});

    list.add_atom(0.0, 0.0, 0.0); // feels the next atom only
    list.add_atom(1.0, 0.0, 0.0); // is pulled equally in both directions
    list.add_atom(2.0, 0.0, 0.0); // feels the previous atom only

    calc_forces_internal(list, TestFF);

    const auto force = -24.0;

    const vector<double> expected {
        force, 0.0, 0.0,
        0.0, 0.0, 0.0,
        -force, 0.0, 0.0
    };

    ASSERT_EQ_VEC(list.fs, expected, "force calculation does not add correctly");
)

ADD_TEST(test_calc_force_between_two_lists,
    const RVec size {0.0, 0.0, 0.0};

    CellList list1 (1, RVec {0.0, 0.0, 0.0}, size);
    list1.add_atom(0.0, 0.0, 0.0);

    CellList list2 (1, RVec {1.0, 0.0, 0.0}, size); // shifted 1 along x
    list2.add_atom(0.0, 0.0, 0.0);

    calc_forces_cell_to_cell(list1, list2, TestFF);

    const auto force = -24.0; // epsilon = 24, sigma = 1, dr = 1 (along x only)

    const vector<double> expected1 {force, 0.0, 0.0};
    const vector<double> expected2 {-force, 0.0, 0.0};

    ASSERT_EQ_VEC(list1.fs, expected1, "forces are not calculated correctly");
    ASSERT_EQ_VEC(list2.fs, expected2, "forces are not calculated correctly");
)

ADD_TEST(test_calc_force_internal_and_between_cells_consistency_check,
    const RVec size {1.0, 1.0, 1.0};

    CellList list1 (1, RVec {0.0, 0.0, 0.0}, size);
    list1.add_atom(0.0, 0.0, 0.0);

    CellList list2 (1, RVec {0.1, 0.2, 0.3}, size);
    list2.add_atom(0.0, 0.0, 0.0);

    calc_forces_cell_to_cell(list1, list2, TestFF);

    CellList expect (2, RVec {0.0, 0.0, 0.0}, size);
    expect.add_atom(0.0, 0.0, 0.0);
    expect.add_atom(0.1, 0.2, 0.3);

    calc_forces_internal(expect, TestFF);

    ASSERT_EQ(expect.fs[0], list1.fs[0], "force calc 1 is not consistent");
    ASSERT_EQ(expect.fs[1], list1.fs[1], "force calc 2 is not consistent");
    ASSERT_EQ(expect.fs[2], list1.fs[2], "force calc 3 is not consistent");
    ASSERT_EQ(expect.fs[3], list2.fs[0], "force calc 4 is not consistent");
    ASSERT_EQ(expect.fs[4], list2.fs[1], "force calc 5 is not consistent");
    ASSERT_EQ(expect.fs[5], list2.fs[2], "force calc 6 is not consistent");
)

ADD_TEST(test_calc_force_for_a_wall_if_outside_box,
    const std::string title = "Test";
    const RVec box_size {1.0, 1.0, 1.0};
    const RVec cell_size {1.0, 1.0, 1.0};

    auto system = System(title, box_size);
    system.cell_size = cell_size;
    system.shape = IVec {1, 1, 1};

    // The list has two atoms outside the box (above/below the box)
    // and one inside which should not be affected:
    // the list is also translated by 0.5 along all directions, which
    // has to be compensated for!
    auto list = CellList(3, RVec {0.5, 0.5, 0.5}, cell_size);
    list.add_atom(-1.5, -2.5, -3.5); // at (-1, -2, -3) in system coords
    list.add_atom(1.5, 2.5, 3.5); // at (2, 3, 4) in system coords
    list.add_atom(0.0, 0.0, 0.0); // at (0.5, 0.5, 0.5) in system coords

    system.cell_lists.push_back(list);

    // Calculate the initial forces (no internal interactions)
    calc_wall_forces(system, TestFF);

    const auto& forces = system.cell_lists[0].fs;

    // Atom 1 should move in positive directions
    ASSERT_EQ(forces[XX], TestFF.wall_constant * 1.0,
        "the force of the first atom was not added correctly");
    ASSERT_EQ(forces[YY], TestFF.wall_constant * 2.0,
        "the force of the first atom was not added correctly");
    ASSERT_EQ(forces[ZZ], TestFF.wall_constant * 3.0,
        "the force of the first atom was not added correctly");

    // Atom 2 should be reversed
    ASSERT_EQ(forces[NDIM + XX], -TestFF.wall_constant * 1.0,
        "the force of the second atom was not added correctly");
    ASSERT_EQ(forces[NDIM + YY], -TestFF.wall_constant * 2.0,
        "the force of the second atom was not added correctly");
    ASSERT_EQ(forces[NDIM + ZZ], -TestFF.wall_constant * 3.0,
        "the force of the second atom was not added correctly");

    // Atom 3 should be unaffected
    ASSERT_EQ(forces[2 * NDIM + XX], 0.0,
        "the third atom was unexpectedly affected by the wall");
    ASSERT_EQ(forces[2 * NDIM + YY], 0.0,
        "the third atom was unexpectedly affected by the wall");
    ASSERT_EQ(forces[2 * NDIM + ZZ], 0.0,
        "the third atom was unexpectedly affected by the wall");
)

ADD_TEST(test_update_velocities,
    CellList list (2, RVec {0.0, 0.0, 0.0}, RVec{1.0, 1.0, 1.0});

    list.add_atom(0.0, 0.0, 0.0);
    list.add_atom(1.0, 0.0, 0.0);

    const vector<real> fs      {0.0, 1.0, 2.0, 3.0,  4.0,  5.0};
    const vector<real> fs_prev {6.0, 7.0, 8.0, 9.0, 10.0, 11.0};

    // Set the current velocity to ensure that the updated velocity
    // takes it into account and does not use it as zero.
    const real old_velocity = 1.0;
    const vector<real> vs (fs.size(), old_velocity);

    // Calculate the expected values
    auto iter = fs.cbegin();
    auto iter_prev = fs_prev.cbegin();
    vector<real> expected;

    while ((iter != fs.cend()) && (iter_prev != fs_prev.cend()))
    {
        // Calculate acceleration as the mean of the current and previous forces
        // according to the Velocity Verlet algorithm
        const auto mean_force = (*iter++ + *iter_prev++) / 2.0;
        const auto acceleration = mean_force / TestFF.mass;

        const auto new_velocity = old_velocity + acceleration * TestOpts.dt;
        expected.push_back(new_velocity);
    }

    list.vs = vs;
    list.fs = fs;
    list.fs_prev = fs_prev;

    update_velocities_cell(list, TestFF, TestOpts);

    ASSERT_EQ_VEC(list.vs, expected, "velocities are not updated correctly");
)

ADD_TEST(test_update_positions,
    CellList list (1, RVec {0.0, 0.0, 0.0}, RVec{1.0, 1.0, 1.0});

    list.add_atom(1.0, 2.0, 3.0);

    constexpr real force = 5e-1;
    constexpr real velocity = 5.0 * force;
    const real expected_dx =
        velocity * TestOpts.dt + 0.5 * (force / TestFF.mass) * TestOpts.dt2;


    vector<real> expected;
    for (const auto& x0 : list.xs)
    {
        expected.push_back(x0 + expected_dx);
    }

    const vector<real> fs (list.xs.size(), force);
    const vector<real> vs (list.xs.size(), velocity);
    list.vs = vs;
    list.fs = fs;

    update_positions_cell(list, TestFF, TestOpts);

    ASSERT_EQ_VEC(list.xs, expected, "positions are not updated correctly");
)

// Move the current forces in list.fs to list.fs_prev and then set
// all values in list.fs to 0.
static void reset_forces_cell(CellList& list)
{
    list.fs.swap(list.fs_prev);
    list.fs.assign(list.fs.size(), 0.0);
}

ADD_TEST(test_velocity_verlet_step_single_list,
    const std::string title = "Test";
    const RVec box_size {1.0, 1.0, 1.0};
    auto system = System(title, box_size);

    auto list = CellList(2, RVec {0.0, 0.0, 0.0}, box_size);
    list.add_atom(0.0, 0.0, 0.0);
    list.add_atom(1.0, 0.0, 0.0);

    // Calculate the initial forces
    calc_forces_internal(list, TestFF);

    // Add an independent copy of the initial state to the system
    system.cell_lists.push_back(list);

    // Run two steps through the Velocity Verlet scheme manually and using
    // the function to ensure that the functionality is correct.

    // Currently: x0, v0, f0, _ (prev forces empty)
    update_positions_cell(list, TestFF, TestOpts);  // x1, v0, f0, _
    reset_forces_cell(list);                        // x1, v0,  _, f0
    calc_forces_internal(list, TestFF);            // x1, v0, f1, f0
    calc_wall_forces_in_list(list, box_size, TestFF.wall_constant);
    update_velocities_cell(list, TestFF, TestOpts); // x1, v1, f1, f0
    update_positions_cell(list, TestFF, TestOpts);  // x2, v1, f1, f0
    reset_forces_cell(list);                        // x2, v1,  _, f1
    calc_forces_internal(list, TestFF);            // x2, v1, f2, f1
    calc_wall_forces_in_list(list, box_size, TestFF.wall_constant);
    update_velocities_cell(list, TestFF, TestOpts); // x2, v2, f2, f1

    // With the function
    Benchmark bench;
    const auto mpi_comm = get_single_rank_mpi_comm(system);

    run_velocity_verlet(system, bench, mpi_comm, TestFF, TestOpts);
    run_velocity_verlet(system, bench, mpi_comm, TestFF, TestOpts);

    ASSERT_EQ_VEC(system.cell_lists[0].xs, list.xs,
        "positions were not updated correctly");
    ASSERT_EQ_VEC(system.cell_lists[0].vs, list.vs,
        "velocities were not updated correctly");
    ASSERT_EQ_VEC(system.cell_lists[0].fs, list.fs,
        "forces were not updated correctly");
    ASSERT_EQ_VEC(system.cell_lists[0].fs_prev, list.fs_prev,
        "previous forces were not updated correctly");
)

ADD_TEST(test_velocity_verlet_step_list_with_a_neighbour,
    const std::string title = "Test";
    const RVec box_size {2.0, 1.0, 1.0};
    const RVec cell_size {1.0, 1.0, 1.0};

    auto system = System(title, box_size);
    system.cell_size = cell_size;
    system.shape = IVec {2, 1, 1};

    // Two lists separated by 1 along the x axis
    auto list1 = CellList(2, RVec {0.0, 0.0, 0.0}, cell_size);
    auto list2 = CellList(2, RVec {1.0, 0.0, 0.0}, cell_size);
    list1.add_atom(0.6, 0.0, 0.0);
    list2.add_atom(0.4, 0.0, 0.0);

    // Calculate the initial forces (no internal interactions)
    calc_forces_cell_to_cell(list1, list2, TestFF);

    // Copy them to the system as index 0 and 1
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);

    // list1 will interact with the second list
    system.cell_lists[0].to_neighbours.push_back(1);

    // Run through one step of the Velocity Verlet scheme manually
    // for both lists
    update_positions_cell(list1, TestFF, TestOpts);
    update_positions_cell(list2, TestFF, TestOpts);
    reset_forces_cell(list1);
    reset_forces_cell(list2);
    calc_forces_cell_to_cell(list1, list2, TestFF); // still the only interaction
    update_velocities_cell(list1, TestFF, TestOpts);
    update_velocities_cell(list2, TestFF, TestOpts);

    // With the function
    Benchmark bench;
    const auto mpi_comm = get_single_rank_mpi_comm(system);

    run_velocity_verlet(system, bench, mpi_comm, TestFF, TestOpts);

    ASSERT_EQ_VEC(system.cell_lists[0].xs, list1.xs,
        "positions were not updated correctly in list1");
    ASSERT_EQ_VEC(system.cell_lists[1].xs, list2.xs,
        "positions were not updated correctly in list2");
    ASSERT_EQ_VEC(system.cell_lists[0].vs, list1.vs,
        "velocities were not updated correctly in list1");
    ASSERT_EQ_VEC(system.cell_lists[1].vs, list2.vs,
        "velocities were not updated correctly in list2");
    ASSERT_EQ_VEC(system.cell_lists[0].fs, list1.fs,
        "forces were not updated correctly in list1");
    ASSERT_EQ_VEC(system.cell_lists[1].fs, list2.fs,
        "forces were not updated correctly in list2");
    ASSERT_EQ_VEC(system.cell_lists[0].fs_prev, list1.fs_prev,
        "previous forces were not updated correctly in list1");
    ASSERT_EQ_VEC(system.cell_lists[1].fs_prev, list2.fs_prev,
        "previous forces were not updated correctly in list2");
)

ADD_TEST(test_velocity_verlet_step_list_with_no_neighbours_does_nothing,
    const std::string title = "Test";
    const RVec box_size {2.0, 1.0, 1.0};
    auto system = System(title, box_size);

    system.cell_size = RVec {1.0, 1.0, 1.0};
    system.shape = IVec {2, 1, 1};

    // Two lists separated by 1 along the x axis
    auto list1 = CellList(2, RVec {0.0, 0.0, 0.0}, system.cell_size);
    auto list2 = CellList(2, RVec {1.0, 0.0, 0.0}, system.cell_size);
    list1.add_atom(0.0, 0.0, 0.0);
    list2.add_atom(0.0, 0.0, 0.0);

    // Copy them to the system as index 0 and 1
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);

    Benchmark bench;
    const auto mpi_comm = get_single_rank_mpi_comm(system);

    run_velocity_verlet(system, bench, mpi_comm, TestFF, TestOpts);

    const vector<real> zeroes (list1.xs.size(), 0.0);

    ASSERT_EQ_VEC(system.cell_lists[0].xs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
    ASSERT_EQ_VEC(system.cell_lists[0].vs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
    ASSERT_EQ_VEC(system.cell_lists[0].fs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
    ASSERT_EQ_VEC(system.cell_lists[1].xs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
    ASSERT_EQ_VEC(system.cell_lists[1].vs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
    ASSERT_EQ_VEC(system.cell_lists[1].fs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
)

ADD_TEST(test_velocity_verlet_for_system_includes_wall_force_calc,
    const std::string title = "Test";
    const RVec box_size {1.0, 1.0, 1.0};
    const RVec cell_size {1.0, 1.0, 1.0};

    auto system = System(title, box_size);
    system.cell_size = cell_size;
    system.shape = IVec {1, 1, 1};

    // The system has a single atom outside of it, for which a wall
    // force will be added!
    auto list = CellList(3, RVec {0.0, 0.0, 0.0}, cell_size);
    list.add_atom(-1.0, -2.0, -3.0);

    system.cell_lists.push_back(list);

    // Calculate the initial forces (no internal interactions)
    Benchmark bench;
    const auto mpi_comm = get_single_rank_mpi_comm(system);

    run_velocity_verlet(system, bench, mpi_comm, TestFF, TestOpts);

    const auto& forces = system.cell_lists[0].fs;

    ASSERT_EQ(forces[XX], TestFF.wall_constant * 1.0,
        "the force of the first atom was not added correctly");
    ASSERT_EQ(forces[YY], TestFF.wall_constant * 2.0,
        "the force of the first atom was not added correctly");
    ASSERT_EQ(forces[ZZ], TestFF.wall_constant * 3.0,
        "the force of the first atom was not added correctly");
)

RUN_TESTS(
    test_calc_force();
    test_calc_force_with_no_atoms_in_cell_works();
    test_calc_force_outside_of_rcut_is_zero();
    test_calc_force_adds_total_force();
    test_calc_force_between_two_lists();
    test_calc_force_for_a_wall_if_outside_box();
    test_calc_force_internal_and_between_cells_consistency_check();
    test_update_velocities();
    test_update_positions();
    test_velocity_verlet_step_single_list();
    test_velocity_verlet_step_list_with_a_neighbour();
    test_velocity_verlet_step_list_with_no_neighbours_does_nothing();
    test_velocity_verlet_for_system_includes_wall_force_calc();
)
