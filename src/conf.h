#include <array>
#include <cstdint>
#include <string>
#include <vector>

#ifndef SYSTEM_CONF_H
#define SYSTEM_CONF_H

enum Direction {
    XX,
    YY,
    ZZ,
    NDIM
};

using real = double;
using RVec = std::array<real, NDIM>;
using IVec = std::array<uint64_t, NDIM>;

// A cubic box of atoms as a part of the whole system.
class Box {
public:
    // Allocate memory for the box and return an empty configuration.
    Box(const uint64_t capacity, const RVec origin, const RVec size);

    /************
    * Functions *
    *************/
    // Return the number of atoms in the box.
    uint64_t num_atoms(void) const { return natoms; };

    // Add an atom with input position to the system.
    void add_atom(const real x, const real y, const real z);

    /************
    * Variables *
    *************/
    std::vector<real> xs; // Positions (x), 1 elem per dimension: relative to the origin
    std::vector<real> vs; // Velocities (v)
    std::vector<real> fs; // Forces (f)
    std::vector<real> fs_prev; // Forces at previous step for Velocity Verlet

    // Box origin.
    RVec origin;

    // Box size.
    RVec size;

    // Neighbouring boxes indices in a collection (ie. a `System` struct)
    // which this box will interact with and add forces *to*. Since the force
    // calculation is symmetric we only have to calculate it once, going
    // *from* a box *to* another. Thus all boxes in this collection should
    // not have this box's index in its corresponding `to_neighbours` object.
    std::vector<size_t> to_neighbours;

private:
    uint64_t natoms;
};

class System {
public:
    // Prepare a container for the system.
    System(const std::string& title, const RVec box_size);

    /************
    * Functions *
    *************/
    // Return the number of atoms in the system.
    uint64_t num_atoms(void) const;

    /************
    * Variables *
    *************/
    // The system is split into (maybe) several boxes containing the atoms.
    std::vector<Box> boxes;

    // System box size.
    RVec box_size;

    // System shape.
    IVec shape;

    // Title of system.
    std::string title;
};

// Read a configuration from a Gromos formatted file.
System read_conf_from_grofile(const std::string& filename);

// Write the configuration to a Gromos formatted file.
void write_conf_to_grofile(const System& system, const std::string& path);

#endif // SYSTEM_CONF_H
