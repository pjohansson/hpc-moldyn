#include <fstream>
#include <iomanip>

#include "conf.h"

System::System(const std::string& title, const RVec box_size)
    :box_size { box_size },
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
    system.boxes.push_back(box);

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
