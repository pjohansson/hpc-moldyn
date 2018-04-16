#include "tests/utils.h"

#include "src/analytics.cpp"
#include "src/conf.h"
#include "src/integrator.cpp"
#include "src/params.h"

using namespace std;

constexpr ForceField TestFF (
    1.0, // epsilon
    1.0, // sigma
    10.0, // rcut
    2.0  // mass
);

System get_system()
{
    CellList list (1, RVec {0.0, 0.0, 0.0}, RVec {2.0, 2.0, 2.0});

    list.add_atom(0.0, 0.0, 0.0);
    list.add_atom(0.9, 1.0, 1.1);
    list.add_atom(1.9, 1.9, 1.9);

    const vector<double> velocities = {
        0.0, 1.0, 2.0,
        3.0, 4.0, 5.0,
        6.0, 7.0, 8.0
    };
    list.vs = velocities;

    const std::string title {"title of system"};
    const RVec box_size {2.0, 2.0, 2.0};

    System system {title, box_size};
    system.cell_lists.push_back(list);

    return system;
}

ADD_TEST(test_calc_temperature_using_equipartition_theorem_in_single_cell,
    const auto system = get_system();

    const auto Ekin =
          0.5 * (0.0 * 0.0 + 1.0 * 1.0 + 2.0 * 2.0)
        + 0.5 * (3.0 * 3.0 + 4.0 * 4.0 + 5.0 * 5.0)
        + 0.5 * (6.0 * 6.0 + 7.0 * 7.0 + 8.0 * 8.0);
    const double N_dof = NDIM * 2.0; // One particle (free from some relative particle)
    const auto temp = 2 * Ekin / N_dof;

    Energetics energy;

    ASSERT_EQ(energy.potential.size(), 0,
        "Energetics object incorrectly initialized");
    ASSERT_EQ(energy.temperature.size(), 0,
        "Energetics object incorrectly initialized");

    calculate_system_energetics(energy, system, TestFF);

    ASSERT_EQ(energy.temperature.back(), temp,
        "temperature not calculated correctly");
)

ADD_TEST(test_calc_potential_energy_in_single_cell_system,
    const auto system = get_system();
    const auto& xs = system.cell_lists[0].xs;

    const auto dr01 = calc_distance(xs, 0, 1);
    const auto dr02 = calc_distance(xs, 0, 2);
    const auto dr12 = calc_distance(xs, 1, 2);

    // Just add up the potential
    const double potential = 4.0 * (
        1.0 / std::pow(dr01[NDIM], 6) - 1.0 / std::pow(dr01[NDIM], 3)
        + 1.0 / std::pow(dr02[NDIM], 6) - 1.0 / std::pow(dr02[NDIM], 3)
        + 1.0 / std::pow(dr12[NDIM], 6) - 1.0 / std::pow(dr12[NDIM], 3)
    );

    Energetics energy;
    calculate_system_energetics(energy, system, TestFF);

    ASSERT_EQ(energy.potential.back(), potential,
        "potential not calculated correctly");
)

ADD_TEST(test_calc_energetics_in_multi_cell_system,
    auto system = get_system();
    const auto num_atoms = system.cell_lists[0].num_atoms();

    // Calculate the reference values from the unsplit system
    Energetics energy;
    calculate_system_energetics(energy, system, TestFF);

    // Split the system at 1.0
    create_cell_lists(system, 0.99);

    ASSERT(system.cell_lists.size() > 1, "system was not split correctly");
    ASSERT(static_cast<bool>(system.cell_lists[0].num_atoms() < num_atoms),
        "the first cell still contains all the atoms, which should not happen");

    // Compare to the values afterwards
    calculate_system_energetics(energy, system, TestFF);
    ASSERT_EQ(energy.temperature.back(), energy.temperature.front(),
        "temperature not calculated correctly in divided system");
    ASSERT_EQ(energy.potential.back(), energy.potential.front(),
        "potential energy not calculated correctly in divided system");
)

RUN_TESTS(
    test_calc_temperature_using_equipartition_theorem_in_single_cell();
    test_calc_potential_energy_in_single_cell_system();
    test_calc_energetics_in_multi_cell_system();
)
