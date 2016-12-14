#include "conf.h"

constexpr int dim = 3;

system_conf::system_conf(int capacity)
    :num_atoms{0}
{
    xs.reserve(dim*capacity);
    vs.reserve(dim*capacity);
    fs.reserve(dim*capacity);
}

int system_conf::add_atom(double x, double y, double z)
{
    xs.push_back(x);
    xs.push_back(y);
    xs.push_back(z);

    for (int i = 0; i < dim; ++i) {
        vs.push_back(0.0);
        fs.push_back(0.0);
    }

    ++num_atoms;

    return num_atoms;
}
