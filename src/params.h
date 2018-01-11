#include <cmath>

#include "conf.h"

#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

struct ForceField {
    constexpr ForceField(const real epsilon,
                         const real sigma,
                         const real rcut,
                         const real mass)
    :epsilon { epsilon },
     sigma { sigma },
     mass { mass },
     c6 { 24.0 * epsilon * std::pow(sigma, 6) },
     c12 { 48.0 * epsilon * std::pow(sigma, 12) },
     rcut { rcut },
     rcut2 { std::pow(rcut, 2) }
    {}

    real epsilon,
         sigma,
         mass,
         c6,
         c12,
         rcut,
         rcut2;
};

struct Options {
    constexpr Options(const double dt)
    :dt { dt },
     dt2 { dt * dt }
    {}

    double dt, dt2;
};

constexpr ForceField DefaultFF = ForceField(4.39e-1, 0.28, 1.1, 12e-5);
constexpr Options DefaultOpts = Options(5e-12);

#endif // SIMULATION_PARAMETERS_H
