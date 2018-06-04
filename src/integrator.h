#include "analytics.h"
#include "conf.h"
#include "mpi_impl.h"
#include "params.h"

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

// Run through a simulation step using the Velocity Verlet scheme.
void run_velocity_verlet(System& system,
                         Benchmark& bench,
                         const MPIRank& mpi_comm,
                         const ForceField& ff,
                         const Options& opts);

#endif // INTEGRATOR_H
