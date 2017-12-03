#include <cmath>

#include "conf.h"

#ifndef FORCE_FIELD_H
#define FORCE_FIELD_H

struct ForceField {
    constexpr ForceField(const real epsilon,
                         const real sigma,
                         const real rcut)
    :epsilon { epsilon },
     sigma { sigma },
     c6 { 24.0 * epsilon * std::pow(sigma, 6) },
     c12 { 48.0 * epsilon * std::pow(sigma, 12) },
     rcut { rcut },
     rcut2 { std::pow(rcut, 2) }
    {}

    real epsilon,
         sigma,
         c6,
         c12,
         rcut,
         rcut2;
}

constexpr DefaultFF = ForceField(1.0, 1.0, 1.0);

#endif // FORCE_FIELD_H
