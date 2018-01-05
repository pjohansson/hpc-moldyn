#include <vector>

#include "tests/utils.h"

#include "src/conf.h"
#include "src/params.h"
#include "src/integrator.cpp"

using namespace std;

constexpr ForceField TestFF (
    1.0, // epsilon
    1.0, // sigma
    1.1, // rcut
    2.0  // mass
);

constexpr Options TestOpts { 5e-3 };

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

ADD_TEST(test_update_velocities,
    Box box (2, RVec {0.0, 0.0, 0.0}, RVec{1.0, 1.0, 1.0});

    box.add_atom(0.0, 0.0, 0.0);
    box.add_atom(1.0, 0.0, 0.0);

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

    box.vs = vs;
    box.fs = fs;
    box.fs_prev = fs_prev;

    update_velocities_box(box, TestFF, TestOpts);

    ASSERT_EQ_VEC(box.vs, expected, "velocities are not updated correctly");
)

ADD_TEST(test_update_positions,
    Box box (1, RVec {0.0, 0.0, 0.0}, RVec{1.0, 1.0, 1.0});

    box.add_atom(1.0, 2.0, 3.0);

    constexpr real force = 5e-1;
    constexpr real velocity = 5.0 * force;
    constexpr real expected_dx =
        velocity * TestOpts.dt + 0.5 * (force / TestFF.mass) * TestOpts.dt2;


    vector<real> expected;
    for (const auto& x0 : box.xs)
    {
        expected.push_back(x0 + expected_dx);
    }

    const vector<real> fs (box.xs.size(), force);
    const vector<real> vs (box.xs.size(), velocity);
    box.vs = vs;
    box.fs = fs;

    update_positions_box(box, TestFF, TestOpts);

    ASSERT_EQ_VEC(box.xs, expected, "positions are not updated correctly");
)

ADD_TEST(test_reset_forces_sets_current_to_previous_and_resets,
    auto box = Box(2, RVec {0.0, 0.0, 0.0}, RVec {1.0, 1.0, 1.0});
    box.add_atom(0.0, 0.0, 0.0);
    box.add_atom(1.0, 0.0, 0.0);

    constexpr real force = 1.0;
    const vector<real> fs (box.xs.size(), force);

    // Also set the previous forces to non-zero values which should be discarded
    constexpr real force_prev = 0.5 * force;
    const vector<real> fs_prev (box.xs.size(), force_prev);

    box.fs = fs;
    box.fs_prev = fs_prev;

    const auto p0 = box.fs.data();
    reset_forces_box(box);
    const auto p1 = box.fs_prev.data();

    const vector<real> zeroes (box.xs.size(), 0.0);

    ASSERT_EQ_VEC(box.fs, zeroes, "forces were not set to zero");
    ASSERT_EQ_VEC(box.fs_prev, fs, "previous forces were not set to the current");
    ASSERT_EQ(p0, p1, "the forces were copied instead of moved");
)

ADD_TEST(test_velocity_verlet_step_single_box,
    const std::string title = "Test";
    const RVec box_size {1.0, 1.0, 1.0};
    auto system = System(title, box_size);

    auto box = Box(2, RVec {0.0, 0.0, 0.0}, box_size);
    box.add_atom(0.0, 0.0, 0.0);
    box.add_atom(1.0, 0.0, 0.0);

    // Calculate the initial forces
    calc_forces_internal(box, TestFF);

    // Add an independent copy of the initial state to the system
    system.boxes.push_back(box);

    // Run two steps through the Velocity Verlet scheme manually and using
    // the function to ensure that the functionality is correct.

    // Currently: x0, v0, f0, _ (prev forces empty)
    update_positions_box(box, TestFF, TestOpts);  // x1, v0, f0, _
    reset_forces_box(box);                        // x1, v0,  _, f0
    calc_forces_internal(box, TestFF);            // x1, v0, f1, f0
    update_velocities_box(box, TestFF, TestOpts); // x1, v1, f1, f0
    update_positions_box(box, TestFF, TestOpts);  // x2, v1, f1, f0
    reset_forces_box(box);                        // x2, v1,  _, f1
    calc_forces_internal(box, TestFF);            // x2, v1, f2, f1
    update_velocities_box(box, TestFF, TestOpts); // x2, v2, f2, f1

    // With the function
    run_velocity_verlet(system, TestFF, TestOpts);
    run_velocity_verlet(system, TestFF, TestOpts);

    ASSERT_EQ_VEC(system.boxes[0].xs, box.xs,
        "positions were not updated correctly");
    ASSERT_EQ_VEC(system.boxes[0].vs, box.vs,
        "velocities were not updated correctly");
    ASSERT_EQ_VEC(system.boxes[0].fs, box.fs,
        "forces were not updated correctly");
    ASSERT_EQ_VEC(system.boxes[0].fs_prev, box.fs_prev,
        "previous forces were not updated correctly");
)

ADD_TEST(test_velocity_verlet_step_box_with_a_neighbour,
    const std::string title = "Test";
    const RVec box_size {2.0, 1.0, 1.0};
    auto system = System(title, box_size);

    // Two boxes separated by 1 along the x axis
    auto box1 = Box(2, RVec {0.0, 0.0, 0.0}, RVec {1.0, 1.0, 1.0});
    auto box2 = Box(2, RVec {1.0, 0.0, 0.0}, RVec {1.0, 1.0, 1.0});
    box1.add_atom(0.0, 0.0, 0.0);
    box2.add_atom(0.0, 0.0, 0.0);

    // Calculate the initial forces (no internal interactions)
    calc_forces_from_to_box(box1, box2, TestFF);

    // Copy them to the system as index 0 and 1
    system.boxes.push_back(box1);
    system.boxes.push_back(box2);

    // box1 will interact with the second box
    system.boxes[0].to_neighbours.push_back(1);

    // Run through one step of the Velocity Verlet scheme manually
    // for both boxes
    update_positions_box(box1, TestFF, TestOpts);
    update_positions_box(box2, TestFF, TestOpts);
    reset_forces_box(box1);
    reset_forces_box(box2);
    calc_forces_from_to_box(box1, box2, TestFF); // still the only interaction
    update_velocities_box(box1, TestFF, TestOpts);
    update_velocities_box(box2, TestFF, TestOpts);

    // With the function
    run_velocity_verlet(system, TestFF, TestOpts);

    ASSERT_EQ_VEC(system.boxes[0].xs, box1.xs,
        "positions were not updated correctly in box1");
    ASSERT_EQ_VEC(system.boxes[1].xs, box2.xs,
        "positions were not updated correctly in box2");
    ASSERT_EQ_VEC(system.boxes[0].vs, box1.vs,
        "velocities were not updated correctly in box1");
    ASSERT_EQ_VEC(system.boxes[1].vs, box2.vs,
        "velocities were not updated correctly in box2");
    ASSERT_EQ_VEC(system.boxes[0].fs, box1.fs,
        "forces were not updated correctly in box1");
    ASSERT_EQ_VEC(system.boxes[1].fs, box2.fs,
        "forces were not updated correctly in box2");
    ASSERT_EQ_VEC(system.boxes[0].fs_prev, box1.fs_prev,
        "previous forces were not updated correctly in box1");
    ASSERT_EQ_VEC(system.boxes[1].fs_prev, box2.fs_prev,
        "previous forces were not updated correctly in box2");
)

ADD_TEST(test_velocity_verlet_step_box_with_no_neighbours_does_nothing,
    const std::string title = "Test";
    const RVec box_size {2.0, 1.0, 1.0};
    auto system = System(title, box_size);

    // Two boxes separated by 1 along the x axis
    auto box1 = Box(2, RVec {0.0, 0.0, 0.0}, RVec {1.0, 1.0, 1.0});
    auto box2 = Box(2, RVec {1.0, 0.0, 0.0}, RVec {1.0, 1.0, 1.0});
    box1.add_atom(0.0, 0.0, 0.0);
    box2.add_atom(0.0, 0.0, 0.0);

    // Copy them to the system as index 0 and 1
    system.boxes.push_back(box1);
    system.boxes.push_back(box2);

    run_velocity_verlet(system, TestFF, TestOpts);

    const vector<real> zeroes (box1.xs.size(), 0.0);

    ASSERT_EQ_VEC(system.boxes[0].xs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
    ASSERT_EQ_VEC(system.boxes[0].vs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
    ASSERT_EQ_VEC(system.boxes[0].fs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
    ASSERT_EQ_VEC(system.boxes[1].xs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
    ASSERT_EQ_VEC(system.boxes[1].vs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
    ASSERT_EQ_VEC(system.boxes[1].fs, zeroes,
        "the neighbours interacted even though they were not set as neighbours");
)

RUN_TESTS(
    test_calc_force();
    test_calc_force_outside_of_rcut_is_zero();
    test_calc_force_adds_total_force();
    test_calc_force_between_two_boxes();
    test_update_velocities();
    test_update_positions();
    test_velocity_verlet_step_single_box();
    test_velocity_verlet_step_box_with_a_neighbour();
    test_velocity_verlet_step_box_with_no_neighbours_does_nothing();
    test_reset_forces_sets_current_to_previous_and_resets();
)
