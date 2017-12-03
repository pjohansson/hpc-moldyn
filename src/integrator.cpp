#include <cmath>

#include "conf.h"

constexpr double epsilon = 1;
constexpr double sigma = 1;
constexpr double sigma6 = std::pow(sigma, 6);
constexpr double sigma12 = std::pow(sigma, 12);

void calc_forces(SystemConf& conf)
{
    for (unsigned i = 0; i < conf.num_atoms(); ++i)
    {
        for (unsigned j = i + 1; j < conf.num_atoms(); ++j)
        {
            RVec dr {0.0, 0.0, 0.0};
            auto dr2 = 0.0;

            for (int k = 0; k < NDIM; ++k)
            {
                const auto dx = conf.xs.at(j * NDIM + k)
                    - conf.xs.at(i * NDIM + k);

                dr[k] = dx;
                dr2 += dx * dx;
            }

            if (dr2 > 0.0)
            {
                for (int k = 0; k < NDIM; ++k)
                {
                    const auto dx = dr[k];

                    if (dx > 0.0)
                    {
                        const auto dx7 = std::pow(dx, 7);
                        const auto dx13 = std::pow(dx, 13);

                        const auto force = 24 * epsilon
                            * (sigma12 / dx13 - sigma6 / dx7);

                        conf.fs[i * NDIM + k] = force;
                        conf.fs[j * NDIM + k] = -force;
                    }
                }
            }
        }
    }
}
