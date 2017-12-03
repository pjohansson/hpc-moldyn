#include "conf.h"
#include "forcefield.h"

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

void calc_forces_internal(Box& box, const ForceField &ff);
void calc_forces_from_to_box(Box& from_box, Box& to_box, const ForceField &ff);

#endif // INTEGRATOR_H
