#include <cmath>

#include "integrator.h"

// Equal to `RVec` but with the distance squared added as a final element.
using dRVec = std::array<real, NDIM + 1>;

static RVec calc_shift_between_boxes(const Box& from_box, const Box& to_box)
{
    RVec shift {0.0, 0.0, 0.0};

    for (int k = 0; k < NDIM; ++k)
    {
        shift[k] = to_box.origin[k] - from_box.origin[k];
    }

    return shift;
}

// Calculate the distance between to atoms of input indices in a box.
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
// boxes. This is a separate function from the internal calculation since
// they should not be confused and I do not like overloading.
static dRVec calc_distance_different_boxes(const std::vector<real>& from_xs,
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
        const auto dr6 = dr2 * dr2 * dr2;
        const auto dr12 = dr6 * dr6;

        const auto force = (ff.c12 / dr12 - ff.c6 / dr6) / dr2;

        for (int k = 0; k < NDIM; ++k)
        {
            force_vec[k] = force * dr[k];
        }
    }

    return force_vec;
}

// Add the forces from internal interactions within a single box.
void calc_forces_internal(Box& box, const ForceField& ff)
{
    for (unsigned i = 0; i < box.num_atoms() - 1; ++i)
    {
        for (unsigned j = i + 1; j < box.num_atoms(); ++j)
        {
            const auto dr = calc_distance(box.xs, i, j);
            const auto force = calc_force_between_atoms(dr, ff);

            for (int k = 0; k < NDIM; ++k)
            {
                box.fs.at(i * NDIM + k) += force[k];
                box.fs.at(j * NDIM + k) -= force[k];
            }
        }
    }
}

// Add the forces from interactions between all atoms from one box to another.
// This is separate from the calculation in an internal box since I do not
// like overloading, and the functions do slightly different things: for
// internal calculations, the interaction matrix for indices is symmetric
// and we avoid duplicate calculations in a suitable manner (atoms only
// interact once with each other), making the calculation N log N. Since
//  box-to-box calculations contain no duplicate atoms, we must always count
// the full N * M interactions.
void calc_forces_from_to_box(Box& from_box, Box& to_box, const ForceField& ff)
{
    const auto shift = calc_shift_between_boxes(from_box, to_box);

    for (unsigned i = 0; i < from_box.num_atoms(); ++i)
    {
        for (unsigned j = 0; j < to_box.num_atoms(); ++j)
        {
            const auto dr = calc_distance_different_boxes(
                from_box.xs, i,
                to_box.xs, j,
                shift
            );

            const auto force = calc_force_between_atoms(dr, ff);

            for (int k = 0; k < NDIM; ++k)
            {
                from_box.fs.at(i * NDIM + k) += force[k];
                to_box.fs.at(j * NDIM + k) -= force[k];
            }
        }
    }
}

// Update all the positions inside a box using the Velocity Verlet
// integration scheme.
void update_positions_box(Box& box, const ForceField &ff, const Options& opts)
{
    auto iter_vs = box.vs.cbegin();
    auto iter_fs = box.fs.cbegin();

    const real divisor = 2.0 * ff.mass;

    for (auto& x : box.xs)
    {
        x += *iter_vs++ * opts.dt + *iter_fs++ * opts.dt2 / divisor;
    }
}

// Update all the velocities inside a box using the Velocity Verlet
// integration scheme.
void update_velocities_box(Box& box, const ForceField &ff, const Options& opts)
{
    auto iter_fs = box.fs.cbegin();
    auto iter_fs_prev = box.fs_prev.cbegin();

    const real avg_divisor = 2.0 * ff.mass;

    for (auto& v : box.vs)
    {
        const auto a = (*iter_fs++ + *iter_fs_prev++) / avg_divisor;
        v += a * opts.dt;
    }
}
