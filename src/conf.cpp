#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "conf.h"

RVec rvec_add(const RVec r1, const RVec r2)
{
    return RVec {
        r1[XX] + r2[XX],
        r1[YY] + r2[YY],
        r1[ZZ] + r2[ZZ]
    };
}

RVec rvec_sub(const RVec r1, const RVec r2)
{
    return RVec {
        r1[XX] - r2[XX],
        r1[YY] - r2[YY],
        r1[ZZ] - r2[ZZ]
    };
}

System::System(const std::string& title, const RVec box_size)
    :box_size { box_size },
     shape { IVec {1, 1, 1} },
     cell_size { RVec {0.0, 0.0, 0.0} },
     title { title }
{
}

uint64_t System::num_atoms() const
{
    uint64_t num = 0;

    for (const auto& list : cell_lists)
    {
        num += list.num_atoms();
    }

    return num;
}

CellList::CellList(const uint64_t capacity, const RVec origin, const RVec size)
    :origin { origin },
     size { size },
     natoms { 0 }
{
    xs.reserve(NDIM * capacity);
    vs.reserve(NDIM * capacity);
    fs.reserve(NDIM * capacity);
    fs_prev.reserve(NDIM * capacity);
}

void CellList::add_atom(const real x, const real y, const real z)
{
    xs.push_back(x);
    xs.push_back(y);
    xs.push_back(z);

    for (int i = 0; i < NDIM; ++i)
    {
        vs.push_back(0.0);
        fs.push_back(0.0);
        fs_prev.push_back(0.0);
    }

    ++natoms;
}

System read_conf_from_grofile(const std::string& path)
{
    std::ifstream ifs { path, std::ifstream::in };

    constexpr size_t buflen = 256;
    std::string buffer (buflen, ' ');

    getline(ifs, buffer);
    const std::string title = buffer;

    getline(ifs, buffer);
    const auto num_atoms = static_cast<uint64_t>(stoi(buffer));

    auto list = CellList(num_atoms, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0});

    for (unsigned i = 0; i < num_atoms; ++i)
    {
        getline(ifs, buffer);

        const auto x = std::stod(buffer.substr(20, 8));
        const auto y = std::stod(buffer.substr(28, 8));
        const auto z = std::stod(buffer.substr(36, 8));

        list.add_atom(x, y, z);
    }

    getline(ifs, buffer);
    const auto dx = std::stod(buffer.substr(0, 10));
    const auto dy = std::stod(buffer.substr(10, 10));
    const auto dz = std::stod(buffer.substr(20, 10));

    const RVec box_size {dx, dy, dz};
    list.size = box_size;

    auto system = System(title, box_size);
    system.cell_lists.push_back(std::move(list));

    return system;
}

void write_conf_to_grofile(const System& system, const std::string& path)
{
    constexpr char ATOM_NAME[2] = "C";
    constexpr char RESIDUE_NAME[4] = "SOL";

    std::ofstream out { path, std::ofstream::out };

    out << system.title << '\n'
        << system.num_atoms() << '\n';

    out.setf(std::ios::fixed);
    out.precision(3);

    uint64_t n = 0;

    for (const auto &list : system.cell_lists)
    {
        for (unsigned i = 0; i < list.num_atoms(); ++i)
        {
            out << std::setw(5) << std::right << n
                << std::setw(5) << std::left << RESIDUE_NAME
                << std::setw(5) << ATOM_NAME
                << std::setw(5) << n
                << std::setw(8) << list.xs.at(i * NDIM + XX)
                << std::setw(8) << list.xs.at(i * NDIM + YY)
                << std::setw(8) << list.xs.at(i * NDIM + ZZ)
                << '\n';

            ++n;
        }
    }

    out << std::setw(9) << std::right << system.box_size[0] << ' '
        << std::setw(9) << system.box_size[1] << ' '
        << std::setw(9) << system.box_size[2] << '\n';
}

static size_t get_index_within_limits(const RVec x0,
                                      const RVec bin_size,
                                      const IVec num_bins,
                                      const Direction axis)
{
    const auto i = static_cast<int>(floor(x0[axis] / bin_size[axis]));

    // Ensure that they are placed in a list inside the system,
    // with minimum index 0 and maximum n - 1
    return static_cast<size_t>(std::max(0, std::min((int) num_bins[axis] - 1, i)));
}

static size_t get_atom_bin_index(const RVec x0,
                                 const RVec cell_size,
                                 const IVec system_shape)
{
    const auto ix = get_index_within_limits(x0, cell_size, system_shape, XX);
    const auto iy = get_index_within_limits(x0, cell_size, system_shape, YY);
    const auto iz = get_index_within_limits(x0, cell_size, system_shape, ZZ);

    return ix * system_shape[YY] * system_shape[ZZ]
        + iy * system_shape[ZZ]
        + iz;
}

void update_cell_lists(System& system)
{
    std::vector<CellList> new_lists;
    new_lists.reserve(system.cell_lists.size());

    for (const auto& list : system.cell_lists)
    {

        new_lists.push_back(
            CellList(list.num_atoms(), list.origin, system.cell_size)
        );

        new_lists.back().to_neighbours = list.to_neighbours;
    }

    for (unsigned from_index = 0;
         from_index < system.cell_lists.size();
         ++from_index)
    {
        auto& list = system.cell_lists.at(from_index);

        for (unsigned i = 0; i < list.num_atoms(); ++i)
        {
            const auto current = i * NDIM;

            const auto x = RVec {
                list.xs.at(current + XX),
                list.xs.at(current + YY),
                list.xs.at(current + ZZ)
            };
            const auto x0 = rvec_add(list.origin, x);

            const auto to_index = get_atom_bin_index(
                x0, system.cell_size, system.shape);
            auto& to_list = new_lists.at(to_index);

            // If the indices are the same we have the exact coordinates
            // and can transfer them without adjusting between the (identical)
            // origins. Minor save of computational time, more importantly
            // it avoids another floating point rounding error.
            if (from_index == to_index)
            {
                std::copy_n(list.xs.cbegin() + current, NDIM,
                            std::back_inserter(to_list.xs));
            }
            else {
                const auto x1 = rvec_sub(x0, to_list.origin);
                std::copy_n(x1.cbegin(), NDIM,
                            std::back_inserter(to_list.xs));
            }

            std::copy_n(list.vs.cbegin() + current, NDIM,
                        std::back_inserter(to_list.vs));
            std::copy_n(list.fs.cbegin() + current, NDIM,
                        std::back_inserter(to_list.fs_prev));
        }
    }

    for (auto& list : new_lists)
    {
        list.xs.shrink_to_fit();
        list.vs.shrink_to_fit();
        list.fs_prev.shrink_to_fit();

        list.fs.resize(list.xs.size(), 0.0);
        list.fs.shrink_to_fit();

        list.update_num_atoms();
    }

    system.cell_lists = std::move(new_lists);
}

void create_cell_lists(System& system, const real rcut)
{
    const real target_size = 2.0 * rcut;

    const auto nx = std::max((uint64_t) 1, static_cast<uint64_t>(
        floor(system.box_size[XX] / target_size)));
    const auto ny = std::max((uint64_t) 1, static_cast<uint64_t>(
        floor(system.box_size[YY] / target_size)));
    const auto nz = std::max((uint64_t) 1, static_cast<uint64_t>(
        floor(system.box_size[ZZ] / target_size)));

    system.shape = IVec {nx, ny, nz};

    const auto dx = system.box_size[XX] / nx;
    const auto dy = system.box_size[YY] / ny;
    const auto dz = system.box_size[ZZ] / nz;
    const auto cell_size = RVec {dx, dy, dz};
    system.cell_size = cell_size;

    std::vector<CellList> split_lists;
    for (unsigned ix = 0; ix < nx; ++ix)
    {
        for (unsigned iy = 0; iy < ny; ++iy)
        {
            for (unsigned iz = 0; iz < nz; ++iz)
            {
                const auto origin = RVec {ix * dx, iy * dy, iz * dz};
                const auto list = CellList(system.num_atoms(), origin, cell_size);

                split_lists.push_back(std::move(list));
            }
        }
    }

    split_lists.shrink_to_fit();

    if (system.cell_lists.empty())
    {
        system.cell_lists = std::move(split_lists);
    }
    else
    {
        // We are recreating the cell list completely and move all atoms
        // from all cell lists into the first cell list, then update the
        // cell lists to get each atom into their proper list. Strictly
        // speaking this is not necessary (it would be better to move them
        // to their correct cells at once) but we currently only do this once
        // so it does not matter.
        auto& first = split_lists.at(0);

        for (const auto& list : system.cell_lists)
        {
            std::move(list.xs.begin(), list.xs.end(),
                      std::back_inserter(first.xs));
            std::move(list.vs.begin(), list.vs.end(),
                      std::back_inserter(first.vs));
            std::move(list.fs.begin(), list.fs.end(),
                      std::back_inserter(first.fs));
            std::move(list.fs_prev.begin(), list.fs_prev.end(),
                      std::back_inserter(first.fs_prev));
        }

        first.update_num_atoms();
        system.cell_lists = std::move(split_lists);

        update_cell_lists(system);
    }
}
