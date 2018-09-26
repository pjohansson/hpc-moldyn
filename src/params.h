#include <cmath>

#include "conf.h"

#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

constexpr double BOLTZ = 1.38064852e-23; // [J/K]

struct ForceField {
    // Constructor for a force field that was not correctly constructed.
    constexpr ForceField(void)
    :is_valid { false },
     epsilon { -1.0 },
     sigma { -1.0 },
     mass { -1.0 },
     rcut { -1.0 },
     rcut2 { -1.0 },
     wall_constant { -1.0 }
    {}

    // Constructor for a valid force field.
    constexpr ForceField(const real epsilon,
                         const real sigma,
                         const real rcut,
                         const real mass,
                         const real wall_constant)
    :is_valid { true },
     epsilon { epsilon },
     sigma { sigma },
     mass { mass },
     rcut { rcut },
     rcut2 { std::pow(rcut, 2) },
     wall_constant { wall_constant }
    {}

    bool is_valid;
    real epsilon,
         sigma,
         mass,
         rcut,
         rcut2,
         wall_constant;
};

struct Options {
    Options()
    :gen_temp { 0.0 },
     energy_calc { 0 },
     num_steps { 0 },
     traj_stride { 0 },
     gen_velocities { false }
    {
        // For Argon ~ 2 fs
        this->set_dt(1e-3);
    }

    void set_dt(const double dt)
    {
        this->dt = dt;
        this->dt2 = dt * dt;
    }

    double dt, dt2,
           gen_temp;
    uint64_t energy_calc,
             num_steps,
             traj_stride;
    bool gen_velocities;
};

// Argon force field
constexpr ForceField ArgonFF = ForceField(
    1.654e-21, // epsilon [J]
    0.340,  // sigma [nm]
    3.0, // rcut = 3 * sigma (non-dimensional)
    39.948, // molecular mass (u),
    300000.0 // [epsilon] (non-dimensional)
);

// Invalid force field, returned by `read_forcefield_file` if it fails.
constexpr ForceField InvalidFF = ForceField();
constexpr ForceField UninitializedFF = ForceField();

bool read_parameter_file(const std::string path, Options& opts);
ForceField read_forcefield_file(const std::string path);

#endif // SIMULATION_PARAMETERS_H
