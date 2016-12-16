#include <fstream>
#include <iomanip>

#include "conf.h"

using namespace std;

SystemConf::SystemConf(const int capacity)
    :box{0.0, 0.0, 0.0},
     title{""},
     natoms{0}
{
    xs.reserve(NDIM*capacity);
    vs.reserve(NDIM*capacity);
    fs.reserve(NDIM*capacity);
}

int SystemConf::add_atom(const double x, const double y, const double z)
{
    xs.push_back(x);
    xs.push_back(y);
    xs.push_back(z);

    for (int i = 0; i < NDIM; ++i) {
        vs.push_back(0.0);
        fs.push_back(0.0);
    }
    ++natoms;

    return natoms;
}

void SystemConf::set_box(const double x, const double y, const double z)
{
    box[XX] = x;
    box[YY] = y;
    box[ZZ] = z;
}

SystemConf read_conf_from_grofile(const string path)
{
    ifstream ifs { path, ifstream::in };

    constexpr size_t buflen = 256;
    string buffer (buflen, ' ');

    getline(ifs, buffer);
    const string title = buffer;
    getline(ifs, buffer);
    const int num_atoms = stoi(buffer);

    SystemConf conf (num_atoms);
    conf.title = title;

    for (int i = 0; i < num_atoms; ++i) {
        getline(ifs, buffer);
        const auto x = stod(buffer.substr(20, 8));
        const auto y = stod(buffer.substr(28, 8));
        const auto z = stod(buffer.substr(36, 8));
        conf.add_atom(x, y, z);
    }

    getline(ifs, buffer);
    const auto x = stod(buffer.substr(0, 10));
    const auto y = stod(buffer.substr(10, 10));
    const auto z = stod(buffer.substr(20, 10));
    conf.set_box(x, y, z);

    return conf;
}

void write_conf_to_grofile(const SystemConf& conf, const string& path)
{
    ofstream out { path, ofstream::out };

    out << conf.title << '\n'
        << conf.num_atoms() << '\n';

    out.setf(ios::fixed);
    out.precision(3);

    for (int i = 0; i < conf.num_atoms(); ++i) {
        out << setw(5) << right << i << setw(5) << left << RESIDUE_NAME << setw(5) << ATOM_NAME << setw(5) << i
            << setw(8) << conf.xs[i*NDIM] << setw(8) << conf.xs[i*NDIM + 1] << setw(8) << conf.xs[i*NDIM + 2]
            << '\n';
    }

    out << setw(9) << right << conf.box[0] << ' '
        << setw(9) << conf.box[1] << ' '
        << setw(9) << conf.box[2] << '\n';
}
