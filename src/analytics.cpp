#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric> // accumulate

#include "analytics.h"

void describe_system_config(const System& system)
{
    std::cerr << "System title: " << system.title << '\n'
              << "Size:  "
                << system.box_size[XX] << " x "
                << system.box_size[YY] << " x "
                << system.box_size[ZZ] << " (nm^3)\n"
              << "Number of atoms:  " << system.num_atoms() << '\n'
              << "Cell list size:  "
                << system.cell_size[XX] << " x "
                << system.cell_size[YY] << " x "
                << system.cell_size[ZZ] << " (nm^3)\n"
              << "Cell list configuration:  "
                << system.shape[XX] << " x "
                << system.shape[YY] << " x "
                << system.shape[ZZ] << '\n'
              << '\n';

    show_atom_cell_list_distribution(system);
}

void show_atom_cell_list_distribution(const System& system)
{
    constexpr auto max_hist_length = 25;

    uint64_t max_atoms = 0;
    for (const auto& list : system.cell_lists)
    {
        if (list.num_atoms() > max_atoms)
        {
            max_atoms = list.num_atoms();
        }
    }

    unsigned n = 0;
    std::string stars(max_hist_length, ' ');

    std::cerr << "Atom distribution in cell lists:\n";

    for (const auto& list : system.cell_lists)
    {
        const auto frac = list.num_atoms() / static_cast<float>(max_atoms);
        const auto num_chars = static_cast<unsigned>(frac * max_hist_length);

        stars.replace(stars.cbegin(), stars.cend(), max_hist_length, ' ');
        stars.replace(0, num_chars, num_chars, '*');

        std::cerr << std::setw(4) << n << ": "
                  << stars << " (" << list.num_atoms() << ")\n";
        ++n;
    }
}

void Benchmark::finalize(void)
{
    simulation_total = std::chrono::system_clock::now() - total_start;
    rest = simulation_total - cell_list_update - force_update
        - position_update - velocity_update;
}

void Benchmark::start_cell_list_update(void)
{
    cell_list_start = std::chrono::system_clock::now();
}

void Benchmark::start_force_update(void)
{
    force_start = std::chrono::system_clock::now();
}

void Benchmark::start_position_update(void)
{
    position_start = std::chrono::system_clock::now();
}

void Benchmark::start_velocity_update(void)
{
    velocity_start = std::chrono::system_clock::now();
}

void Benchmark::stop_cell_list_update(void)
{
    cell_list_update += std::chrono::system_clock::now() - cell_list_start;
}

void Benchmark::stop_force_update(void)
{
    force_update += std::chrono::system_clock::now() - force_start;
}

void Benchmark::stop_position_update(void)
{
    position_update += std::chrono::system_clock::now() - position_start;
}

void Benchmark::stop_velocity_update(void)
{
    velocity_update += std::chrono::system_clock::now() - velocity_start;
}

namespace bench_print_params {
    constexpr size_t name_length = 30;
    constexpr size_t time_length = 10;
    constexpr size_t perc_length = 8;

    constexpr size_t sep_length = 2;
    constexpr size_t line_length = name_length + time_length + perc_length
        + 2 * sep_length;
}

using namespace std::chrono;
using namespace bench_print_params;

static double calc_percentage(const duration<double> time,
                              const duration<double> total)
{
    return 100.0 * static_cast<double>(time / total);
}

static double calc_seconds(const duration<double> time)
{
    const auto ms = static_cast<double>(
        duration_cast<milliseconds>(time).count()
    );

    return ms / 1000.0;
}

static void print_a_timing(const std::string name,
                           const duration<double> time,
                           const duration<double> total)
{
    const auto sep = std::string(sep_length, ' ');

    std::cerr
        << std::setw(name_length) << std::left << name << std::right
        << std::setprecision(1) << std::fixed
        << sep
        << std::setw(time_length) << calc_seconds(time)
        << sep
        << std::setw(perc_length) << calc_percentage(time, total)
        << '\n';
}

void print_benchmark(const Benchmark& bench)
{
    const auto linebreaker = std::string(line_length, '-') + '\n';

    std::cerr << std::setw(name_length)
        << std::left << "Benchmark Timings" << std::right
        << std::setw(time_length + sep_length) << "Time (s)"
        << std::setw(perc_length + sep_length) << "Frac (%)" << '\n'
        << linebreaker;

    print_a_timing("Cell list update",
                   bench.cell_list_update, bench.simulation_total);
    print_a_timing("Force calculation",
                   bench.force_update, bench.simulation_total);
    print_a_timing("Position update",
                   bench.position_update, bench.simulation_total);
    print_a_timing("Velocity update",
                   bench.velocity_update, bench.simulation_total);
    print_a_timing("Rest",
                   bench.rest, bench.simulation_total);

    std::cerr << linebreaker
        << std::setw(name_length) << std::left << "Total walltime"
        << std::setw(time_length + sep_length)
        << std::right << std::setprecision(1) << std::fixed
        << calc_seconds(bench.simulation_total) << '\n';
}

double calculate_system_temperature(const System& system)
{
    double Ekin_total = 0.0;

    for (const auto& cell : system.cell_lists)
    {
        // Use a simple fold with a lambda expression to square all velocities
        // independently, since we only have one mass.
        const auto Ekin = std::accumulate(
            cell.vs.cbegin(),
            cell.vs.cend(),
            0.0,
            [](const real acc, const real value) {
                return acc + std::pow(value, 2);
            }
        );

        Ekin_total += 0.5 * Ekin;
    }

    const auto degrees_of_freedom = NDIM * (system.num_atoms() - 1);
    const auto temperature = 2.0 * Ekin_total / degrees_of_freedom;

    return temperature;
}

double calc_system_potential_energy(const System& system, const ForceField& ff)
{
    double Epot = 0.0;

    for (const auto& list : system.cell_lists)
    {
        const auto& xs = list.xs;

        for (unsigned i = 0; i < list.num_atoms(); ++i)
        {
            for (unsigned j = i + 1; j < list.num_atoms(); ++j)
            {
                const auto dr2 =
                    std::pow(xs[i * NDIM + XX] - xs[j * NDIM + XX], 2)
                    + std::pow(xs[i * NDIM + YY] - xs[j * NDIM + YY], 2)
                    + std::pow(xs[i * NDIM + ZZ] - xs[j * NDIM + ZZ], 2);

                if (dr2 <= ff.rcut2)
                {
                    const auto dr6inv = 1.0 / std::pow(dr2, 3);
                    const auto dr12inv = std::pow(dr6inv, 2);
                    Epot += 4.0 * (dr12inv - dr6inv);
                }
            }
        }

        for (const auto& to_index : list.to_neighbours)
        {
            const auto& to_list = system.cell_lists.at(to_index);
            const auto& to_xs = to_list.xs;
            const auto dr_box = rvec_sub(to_list.origin, list.origin);

            for (unsigned i = 0; i < list.num_atoms(); ++i)
            {
                for (unsigned j = 0; j < to_list.num_atoms(); ++j)
                {
                    const auto dr2 =
                        std::pow(xs[i * NDIM + XX]
                            - to_xs[j * NDIM + XX]
                            - dr_box[XX], 2)
                        + std::pow(xs[i * NDIM + YY]
                            - to_xs[j * NDIM + YY]
                            - dr_box[YY], 2)
                        + std::pow(xs[i * NDIM + ZZ]
                            - to_xs[j * NDIM + ZZ]
                            - dr_box[ZZ], 2);

                    if (dr2 <= ff.rcut2)
                    {
                        const auto dr6inv = 1.0 / std::pow(dr2, 3);
                        const auto dr12inv = std::pow(dr6inv, 2);

                        Epot += 4.0 * (dr12inv - dr6inv);
                    }
                }
            }
        }
    }

    return Epot;
}

void calculate_system_energetics(Energetics& energy,
                                 const System& system,
                                 const ForceField& ff)
{
    energy.potential.push_back(calc_system_potential_energy(system, ff));
    energy.temperature.push_back(calculate_system_temperature(system));
}
