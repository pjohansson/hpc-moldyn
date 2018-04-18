#include <cmath>

#include "conf.h"

#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

constexpr double BOLTZ = 1.38064852e-23; // [J/K]

struct ForceField {
    constexpr ForceField(const real epsilon,
                         const real sigma,
                         const real rcut,
                         const real mass,
                         const real wall_constant)
    :epsilon { epsilon },
     sigma { sigma },
     mass { mass },
     rcut { rcut },
     rcut2 { std::pow(rcut, 2) },
     wall_constant { wall_constant }
    {}

    real epsilon,
         sigma,
         mass,
         rcut,
         rcut2,
         wall_constant;
};

struct Options {
    constexpr Options(const double dt, const unsigned energy_calc)
    :dt { dt },
     dt2 { dt * dt },
     energy_calc { energy_calc }
    {}

    double dt, dt2;
    unsigned energy_calc;
};

// Argon force field
constexpr ForceField ArgonFF = ForceField(
    1.654e-21, // epsilon [J]
    0.340,  // sigma [nm]
    3.0, // rcut = 3 * sigma (non-dimensional)
    39.948, // molecular mass (u),
    300000.0 // [epsilon] (non-dimensional)
);

// For argon: 1e-3 ~ 2 fs time step
constexpr Options DefaultOpts = Options(1e-3, 5);

#endif // SIMULATION_PARAMETERS_H
