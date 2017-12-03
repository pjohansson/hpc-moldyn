#include <fstream>
#include <iomanip>

#include "conf.h"

SystemConf::SystemConf(const uint64_t capacity)
    :box{0.0, 0.0, 0.0},
     title{""},
     natoms{0}
{
    xs.reserve(NDIM * capacity);
    vs.reserve(NDIM * capacity);
    fs.reserve(NDIM * capacity);
}

void SystemConf::add_atom(const real x, const real y, const real z)
{
    xs.push_back(x);
    xs.push_back(y);
    xs.push_back(z);

    for (int i = 0; i < NDIM; ++i)
    {
        vs.push_back(0.0);
        fs.push_back(0.0);
    }

    ++natoms;
}

void SystemConf::set_box(const real x, const real y, const real z)
{
    box[XX] = x;
    box[YY] = y;
    box[ZZ] = z;
}

SystemConf read_conf_from_grofile(const std::string& path)
{
    std::ifstream ifs { path, std::ifstream::in };

    constexpr size_t buflen = 256;
    std::string buffer (buflen, ' ');

    getline(ifs, buffer);
    const std::string title = buffer;

    getline(ifs, buffer);
    const auto num_atoms = static_cast<uint64_t>(stoi(buffer));

    SystemConf conf (num_atoms);
    conf.title = title;

    for (unsigned i = 0; i < num_atoms; ++i)
    {
        getline(ifs, buffer);
        const auto x = std::stod(buffer.substr(20, 8));
        const auto y = std::stod(buffer.substr(28, 8));
        const auto z = std::stod(buffer.substr(36, 8));
        conf.add_atom(x, y, z);
    }

    getline(ifs, buffer);
    const auto x = std::stod(buffer.substr(0, 10));
    const auto y = std::stod(buffer.substr(10, 10));
    const auto z = std::stod(buffer.substr(20, 10));
    conf.set_box(x, y, z);

    return conf;
}

void write_conf_to_grofile(const SystemConf& conf, const std::string& path)
{
    std::ofstream out { path, std::ofstream::out };

    out << conf.title << '\n'
        << conf.num_atoms() << '\n';

    out.setf(std::ios::fixed);
    out.precision(3);

    for (unsigned i = 0; i < conf.num_atoms(); ++i)
    {
        out << std::setw(5) << std::right << i
            << std::setw(5) << std::left << RESIDUE_NAME
            << std::setw(5) << ATOM_NAME
            << std::setw(5) << i
            << std::setw(8) << conf.xs.at(i * NDIM + XX)
            << std::setw(8) << conf.xs.at(i * NDIM + YY)
            << std::setw(8) << conf.xs.at(i * NDIM + ZZ)
            << '\n';
    }

    out << std::setw(9) << std::right << conf.box[0] << ' '
        << std::setw(9) << conf.box[1] << ' '
        << std::setw(9) << conf.box[2] << '\n';
}
