#include <cmath>
#include <fstream>
#include <iomanip>

#include "conf.h"

System::System(const std::string& title, const RVec box_size)
    :box_size { box_size },
     shape { IVec {1, 1, 1} },
     title { title }
{
}

uint64_t System::num_atoms() const
{
    uint64_t num = 0;

    for (const auto& box : boxes)
    {
        num += box.num_atoms();
    }

    return num;
}

Box::Box(const uint64_t capacity, const RVec origin, const RVec size)
    :origin { origin },
     size { size },
     natoms { 0 }
{
    xs.reserve(NDIM * capacity);
    vs.reserve(NDIM * capacity);
    fs.reserve(NDIM * capacity);
    fs_prev.reserve(NDIM * capacity);
}

void Box::add_atom(const real x, const real y, const real z)
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

Atom Box::get_atom(const size_t index) const
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

System read_conf_from_grofile(const std::string& path)
{
    std::ifstream ifs { path, std::ifstream::in };

    constexpr size_t buflen = 256;
    std::string buffer (buflen, ' ');

    getline(ifs, buffer);
    const std::string title = buffer;

    getline(ifs, buffer);
    const auto num_atoms = static_cast<uint64_t>(stoi(buffer));

    auto box = Box(num_atoms, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0});

    for (unsigned i = 0; i < num_atoms; ++i)
    {
        getline(ifs, buffer);

        const auto x = std::stod(buffer.substr(20, 8));
        const auto y = std::stod(buffer.substr(28, 8));
        const auto z = std::stod(buffer.substr(36, 8));

        box.add_atom(x, y, z);
    }

    getline(ifs, buffer);
    const auto dx = std::stod(buffer.substr(0, 10));
    const auto dy = std::stod(buffer.substr(10, 10));
    const auto dz = std::stod(buffer.substr(20, 10));

    const RVec box_size {dx, dy, dz};
    box.size = box_size;

    auto system = System(title, box_size);
    system.boxes.push_back(std::move(box));

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

    for (const auto &box : system.boxes)
    {
        for (unsigned i = 0; i < box.num_atoms(); ++i)
        {
            out << std::setw(5) << std::right << n
                << std::setw(5) << std::left << RESIDUE_NAME
                << std::setw(5) << ATOM_NAME
                << std::setw(5) << n
                << std::setw(8) << box.xs.at(i * NDIM + XX)
                << std::setw(8) << box.xs.at(i * NDIM + YY)
                << std::setw(8) << box.xs.at(i * NDIM + ZZ)
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

void split_system_into_boxes(System& system, const real rcut)
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
    const auto box_size = RVec {dx, dy, dz};

    std::vector<Box> split_boxes;
    for (unsigned ix = 0; ix < nx; ++ix)
    {
        for (unsigned iy = 0; iy < ny; ++iy)
        {
            for (unsigned iz = 0; iz < nz; ++iz)
            {
                const auto origin = RVec {ix * dx, iy * dy, iz * dz};
                const auto box = Box(system.num_atoms(), origin, box_size);

                split_boxes.push_back(std::move(box));
            }
        }
    }

    for (const auto& box : system.boxes)
    {
        for (unsigned i = 0; i < box.num_atoms(); ++i)
        {
            const auto atom = box.get_atom(i);
            const auto x0 = box.origin[XX] + atom.xs[XX];
            const auto y0 = box.origin[YY] + atom.xs[YY];
            const auto z0 = box.origin[ZZ] + atom.xs[ZZ];

            // Ensure that they are placed in a box inside the system,
            // with minimum index 0 and maximum n - 1
            const auto ix = get_index_within_limits(x0, dx, nx);
            const auto iy = get_index_within_limits(y0, dy, ny);
            const auto iz = get_index_within_limits(z0, dz, nz);

            const auto to_index = ix * ny * nz + iy * nz + iz;

            const auto& origin = split_boxes.at(to_index).origin;
            split_boxes.at(to_index).add_atom(
                x0 - origin[XX], y0 - origin[YY], z0 - origin[ZZ]
            );
        }
    }

    for (auto& box : split_boxes)
    {
        box.xs.shrink_to_fit();
        box.vs.shrink_to_fit();
        box.fs.shrink_to_fit();
        box.fs_prev.shrink_to_fit();
    }

    system.boxes = std::move(split_boxes);
}
