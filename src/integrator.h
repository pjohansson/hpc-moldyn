#include "conf.h"
#include "forcefield.h"

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

void calc_forces_internal(Box& box, const ForceField &ff);
void calc_forces_from_to_box(Box& from_box, Box& to_box, const ForceField &ff);
void update_positions_box(Box& box, const ForceField &ff, const real dt);
void update_velocities_box(Box& box, const ForceField &ff, const real dt);

#endif // INTEGRATOR_H
