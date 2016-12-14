#include <vector>

struct system_conf {
    system_conf(int capacity);

    std::vector<double> xs; // positions (x), 1 elem per dimension
    std::vector<double> vs; // velocities (v)
    std::vector<double> fs; // forces (f)

    int num_atoms;

    int add_atom(double x, double y, double z);
};
