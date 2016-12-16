#include <array>
#include <string>
#include <vector>

#ifndef SYSTEM_CONF_H
#define SYSTEM_CONF_H

enum Direction {
    XX,
    YY,
    ZZ,
    NDIM
};

constexpr char ATOM_NAME[2] = "C";
constexpr char RESIDUE_NAME[4] = "SOL";

class SystemConf {
public:
    SystemConf(int capacity);

    std::vector<double> xs; // positions (x), 1 elem per dimension
    std::vector<double> vs; // velocities (v)
    std::vector<double> fs; // forces (f)

    std::array<double, NDIM> box;
    std::string title;

    int num_atoms() const { return natoms; };

    int add_atom(double x, double y, double z);
    void set_box(double x, double y, double z);

private:
    int natoms;
};

SystemConf read_conf_from_grofile(const std::string filename);
void write_conf_to_grofile(const SystemConf& conf, const std::string& path);

#endif // SYSTEM_CONF_H
