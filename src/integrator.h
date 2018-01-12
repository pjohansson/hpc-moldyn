#include "analytics.h"
#include "conf.h"
#include "params.h"

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

// Run through a simulation step using the Velocity Verlet scheme.
void run_velocity_verlet(System& system,
                         Benchmark& bench,
                         const ForceField& ff,
                         const Options& opts);

#endif // INTEGRATOR_H
