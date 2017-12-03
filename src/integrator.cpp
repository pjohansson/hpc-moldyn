#include <cmath>

#include "conf.h"
#include "forcefield.h"

void calc_forces_internal(SystemConf& conf, const ForceField &ff)
{
    for (unsigned i = 0; i < conf.num_atoms() - 1; ++i)
    {
        for (unsigned j = i + 1; j < conf.num_atoms(); ++j)
        {
            RVec dr_comp {0.0, 0.0, 0.0};
            auto dr2 = 0.0;

            for (int k = 0; k < NDIM; ++k)
            {
                const auto dx = conf.xs.at(j * NDIM + k)
                    - conf.xs.at(i * NDIM + k);

                dr_comp[k] = dx;
                dr2 += dx * dx;
            }

            if (dr2 > 0.0 && dr2 <= ff.rcut2)
            {
                const auto dr6 = dr2 * dr2 * dr2;
                const auto dr12 = dr6 * dr6;

                const auto force = 24 * ff.epsilon
                    * (2.0 * ff.sigma12 / dr12 - ff.sigma6 / dr6);

                for (int k = 0; k < NDIM; ++k)
                {
                    const auto dir_force = force * dr_comp[k] / dr2;

                    conf.fs.at(i * NDIM + k) += dir_force;
                    conf.fs.at(j * NDIM + k) -= dir_force;
                }
            }
        }
    }
}
