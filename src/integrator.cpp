#include <cmath>

#include "integrator.h"

// Equal to `RVec` but with the distance squared added as a final element.
using dRVec = std::array<real, NDIM + 1>;

// Simply return the difference between the two cell lists origins.
static RVec calc_shift_between_cells(const CellList& from_list, const CellList& to_list)
{
    return rvec_sub(to_list.origin, from_list.origin);
}

// Calculate the distance between to atoms of input indices in a cell list.
static dRVec calc_distance(const std::vector<real>& xs,
                           const size_t             from,
                           const size_t             to)
{
    dRVec dr {0.0, 0.0, 0.0, 0.0};

    for (int k = 0; k < NDIM; ++k)
    {
        dr[k] = xs.at(to * NDIM + k) - xs.at(from * NDIM + k);
        dr[NDIM] += dr[k] * dr[k];
    }

    return dr;
}

// Calculate the distance between two atoms of input indices in different
// cells. This is a separate function from the internal calculation since
// they should not be confused and I do not like overloading.
static dRVec calc_distance_different_cells(const std::vector<real>& from_xs,
                                           const size_t             i1,
                                           const std::vector<real>& to_xs,
                                           const size_t             i2,
                                           const RVec&              shift)
{
    dRVec dr {0.0, 0.0, 0.0, 0.0};

    for (int k = 0; k < NDIM; ++k)
    {
        dr[k] = to_xs.at(i2 * NDIM + k) - from_xs.at(i1 * NDIM + k) + shift[k];
        dr[NDIM] += dr[k] * dr[k];
    }

    return dr;
}

static RVec calc_force_between_atoms(const dRVec&      dr,
                                     const ForceField& ff)
{
    RVec force_vec {0.0, 0.0, 0.0};
    const auto dr2 = dr[NDIM];

    if (dr2 > 0.0 && dr2 <= ff.rcut2)
    {
        const auto dr2inv = 1.0 / dr2;
        const auto dr8inv = std::pow(dr2inv, 4);
        const auto dr14inv = std::pow(dr2inv, 7);

        const auto force = 48.0 * (dr14inv - 0.5 * dr8inv);

        for (int k = 0; k < NDIM; ++k)
        {
            force_vec[k] = force * dr[k];
        }
    }

    return force_vec;
}

// Add the forces from internal interactions within a cell list.
static void calc_forces_internal(CellList& list, const ForceField& ff)
{
    if (list.num_atoms() == 0)
    {
        return;
    }

    for (unsigned i = 0; i < list.num_atoms() - 1; ++i)
    {
        for (unsigned j = i + 1; j < list.num_atoms(); ++j)
        {
            const auto dr = calc_distance(list.xs, i, j);
            const auto force = calc_force_between_atoms(dr, ff);

            for (int k = 0; k < NDIM; ++k)
            {
                list.fs.at(i * NDIM + k) += force[k];
                list.fs.at(j * NDIM + k) -= force[k];
            }
        }
    }
}

// Add the forces from interactions between all atoms from one cell list to another.
// This is separate from the calculation in an internal cell list since I do not
// like overloading, and the functions do slightly different things: for
// internal calculations, the interaction matrix for indices is symmetric
// and we avoid duplicate calculations in a suitable manner (atoms only
// interact once with each other), making the calculation N log N. Since
//  list-to-list calculations contain no duplicate atoms, we must always count
// the full N * M interactions.
static void calc_forces_cell_to_cell(CellList& from_list,
                                     CellList& to_list,
                                     const ForceField& ff)
{
    const auto shift = calc_shift_between_cells(from_list, to_list);

    for (unsigned i = 0; i < from_list.num_atoms(); ++i)
    {
        for (unsigned j = 0; j < to_list.num_atoms(); ++j)
        {
            const auto dr = calc_distance_different_cells(
                from_list.xs, i,
                to_list.xs, j,
                shift
            );

            const auto force = calc_force_between_atoms(dr, ff);

            for (int k = 0; k < NDIM; ++k)
            {
                from_list.fs.at(i * NDIM + k) += force[k];
                to_list.fs.at(j * NDIM + k) -= force[k];
            }
        }
    }
}

// Update all the positions inside a cell list using the Velocity Verlet
// integration scheme.
static void update_positions_cell(CellList& list,
                                  const ForceField &ff,
                                  const Options& opts)
{
    auto iter_vs = list.vs.cbegin();
    auto iter_fs = list.fs.cbegin();

    const real divisor = 2.0 * ff.mass;

    for (auto& x : list.xs)
    {
        x += *iter_vs++ * opts.dt + *iter_fs++ * opts.dt2 / divisor;
    }
}

// Update all the velocities inside a cell list using the Velocity Verlet
// integration scheme.
static void update_velocities_cell(CellList& list,
                                   const ForceField &ff,
                                   const Options& opts)
{
    auto iter_fs = list.fs.cbegin();
    auto iter_fs_prev = list.fs_prev.cbegin();

    const real avg_divisor = 2.0 * ff.mass;

    for (auto& v : list.vs)
    {
        const auto a = (*iter_fs++ + *iter_fs_prev++) / avg_divisor;
        v += a * opts.dt;
    }
}

void run_velocity_verlet(System& system,
                         Benchmark& bench,
                         const ForceField& ff,
                         const Options& opts)
{
    bench.start_position_update();
    for (auto& list : system.cell_lists)
    {
        update_positions_cell(list, ff, opts);
    }
    bench.stop_position_update();

    bench.start_cell_list_update();
    update_cell_lists(system);
    bench.stop_cell_list_update();

    // The force calculation is in a separate iteration since
    // it is not local to each cell list, which resetting the forces
    // in each iteration will mess up.
    bench.start_force_update();
    for (auto& list : system.cell_lists)
    {
        calc_forces_internal(list, ff);
        for (const auto& i : list.to_neighbours)
        {
            calc_forces_cell_to_cell(list, system.cell_lists.at(i), ff);
        }
    }
    bench.stop_force_update();

    bench.start_velocity_update();
    for (auto& list : system.cell_lists)
    {
        update_velocities_cell(list, ff, opts);
    }
    bench.stop_velocity_update();
}
