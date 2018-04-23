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
    Options()
    :dt { 1e-3 }, // For argon: 2 fs
     dt2 { dt * dt },
     gen_temp { 0.0 },
     energy_calc { 50 },
     num_steps { 1000 },
     traj_stride { 50 },
     gen_velocities { false },
     verbose { false }
    {}

    void set_dt(const double dt)
    {
        this->dt = dt;
        this->dt2 = dt * dt;
    }

    double dt, dt2,
           gen_temp;
    unsigned energy_calc,
             num_steps,
             traj_stride;
    std::string input_conf,
                output_conf;
    bool gen_velocities,
         verbose;
};

// Argon force field
constexpr ForceField ArgonFF = ForceField(
    1.654e-21, // epsilon [J]
    0.340,  // sigma [nm]
    3.0, // rcut = 3 * sigma (non-dimensional)
    39.948, // molecular mass (u),
    300000.0 // [epsilon] (non-dimensional)
);

bool read_parameter_file(const std::string path, Options& opts);

#endif // SIMULATION_PARAMETERS_H
