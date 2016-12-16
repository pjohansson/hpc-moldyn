#include <array>
#include <cmath>

#include "conf.h"

using namespace std;

constexpr double epsilon = 1;
constexpr double sigma = 1;
constexpr double sigma6 = pow(sigma, 6);
constexpr double sigma12 = pow(sigma, 12);


void calc_forces(SystemConf& conf)
{
    for (int i = 0; i < conf.num_atoms(); ++i) {
        for (int j = i+1; j < conf.num_atoms(); j++) {
            array<double, NDIM> dxs {0.0, 0.0, 0.0};
            auto dx2 = 0.0;
            for (int k = 0; k < NDIM; ++k) {
                auto dx = conf.xs[j*NDIM + k] - conf.xs[i*NDIM + k];
                dxs[k] = dx;
                dx2 += dx*dx;
            }

            if (dx2 > 0.0) {
                for (int k = 0; k < NDIM; ++k) {
                    auto dx = dxs[k];
                    if (dx > 0.0) {
                        auto dx7 = pow(dx, 7);
                        auto dx13 = pow(dx, 13);

                        auto force = 24*epsilon*(sigma12/dx13 - sigma6/dx7);
                        conf.fs[i*NDIM + k] = force;
                        conf.fs[j*NDIM + k] = -force;
                    }
                }
            }
        }
    }
}
