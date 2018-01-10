#include <cmath>
#include <fstream>
#include <iomanip>

#include "conf.h"

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

Atom CellList::get_atom(const size_t index) const
{
    const Atom atom {
        RVec {
            xs.at(NDIM * index + XX),
            xs.at(NDIM * index + YY),
            xs.at(NDIM * index + ZZ)
        },
        RVec {
            vs.at(NDIM * index + XX),
            vs.at(NDIM * index + YY),
            vs.at(NDIM * index + ZZ)
        },
        RVec {
            fs.at(NDIM * index + XX),
            fs.at(NDIM * index + YY),
            fs.at(NDIM * index + ZZ)
        }
    };

    return atom;
}

void CellList::transfer_data_from(const CellList& from_list, const char fields)
{
    natoms = from_list.num_atoms();
    const auto capacity = natoms * NDIM;

    if (fields & POSITION)
    {
        xs = from_list.xs;
    }
    else
    {
        xs.resize(capacity);
        xs.assign(capacity, 0.0);
    }

    if (fields & VELOCITY)
    {
        vs = from_list.vs;
    }
    else
    {
        vs.resize(capacity);
        vs.assign(capacity, 0.0);
    }

    if (fields & FORCE)
    {
        fs = from_list.fs;
        fs_prev = from_list.fs_prev;
    }
    else
    {
        fs.resize(capacity);
        fs.assign(capacity, 0.0);
        fs_prev.resize(capacity);
        fs_prev.assign(capacity, 0.0);
    }
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

static size_t get_index_within_limits(const real x0,
                                      const real bin_size,
                                      const size_t num_bins)
{
    const auto i = static_cast<int>(floor(x0 / bin_size));

    return static_cast<size_t>(std::max(0, std::min((int) num_bins - 1, i)));
}

void update_cell_lists(System& system)
{
    std::vector<CellList> new_lists;
    new_lists.reserve(system.cell_lists.size());

    for (const auto& list : system.cell_lists)
    {
        for (unsigned i = 0; i < list.num_atoms(); ++i)
        {
            const auto atom = list.get_atom(i);
            const auto x0 = list.origin[XX] + atom.xs[XX];
            const auto y0 = list.origin[YY] + atom.xs[YY];
            const auto z0 = list.origin[ZZ] + atom.xs[ZZ];

            // Ensure that they are placed in a list inside the system,
            // with minimum index 0 and maximum n - 1
            const auto ix = get_index_within_limits(
                x0, system.cell_size[XX], system.shape[XX]);
            const auto iy = get_index_within_limits(
                y0, system.cell_size[YY], system.shape[YY]);
            const auto iz = get_index_within_limits(
                z0, system.cell_size[ZZ], system.shape[ZZ]);

            const auto to_index = ix * system.shape[YY] * system.shape[ZZ]
                + iy * system.shape[ZZ] + iz;

            const auto& to_origin = new_lists.at(to_index).origin;
            new_lists.at(to_index).add_atom(
                x0 - to_origin[XX], y0 - to_origin[YY], z0 - to_origin[ZZ]
            );
        }
    }

    for (auto& list : new_lists)
    {
        list.xs.shrink_to_fit();
        list.vs.shrink_to_fit();
        list.fs.shrink_to_fit();
        list.fs_prev.shrink_to_fit();
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

    split_lists.at(0).xs = system.cell_lists.at(0).xs;
    split_lists.at(0).vs = system.cell_lists.at(0).vs;
    split_lists.at(0).fs = system.cell_lists.at(0).fs;
    split_lists.at(0).fs_prev = system.cell_lists.at(0).fs_prev;
    split_lists.shrink_to_fit();
    system.cell_lists = std::move(split_lists);

    update_cell_lists(system);
    // for (const auto& list : system.cell_lists)
    // {
    //     for (unsigned i = 0; i < list.num_atoms(); ++i)
    //     {
    //         const auto atom = list.get_atom(i);
    //         const auto x0 = list.origin[XX] + atom.xs[XX];
    //         const auto y0 = list.origin[YY] + atom.xs[YY];
    //         const auto z0 = list.origin[ZZ] + atom.xs[ZZ];
    //
    //         // Ensure that they are placed in a list inside the system,
    //         // with minimum index 0 and maximum n - 1
    //         const auto ix = get_index_within_limits(x0, dx, nx);
    //         const auto iy = get_index_within_limits(y0, dy, ny);
    //         const auto iz = get_index_within_limits(z0, dz, nz);
    //
    //         const auto to_index = ix * ny * nz + iy * nz + iz;
    //
    //         const auto& origin = split_lists.at(to_index).origin;
    //         split_lists.at(to_index).add_atom(
    //             x0 - origin[XX], y0 - origin[YY], z0 - origin[ZZ]
    //         );
    //     }
    // }
    //
    // for (auto& list : split_lists)
    // {
    //     list.xs.shrink_to_fit();
    //     list.vs.shrink_to_fit();
    //     list.fs.shrink_to_fit();
    //     list.fs_prev.shrink_to_fit();
    // }
    //
}
