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

// An atom's position, velocity and force.
struct Atom {
    RVec xs, vs, fs;
};

// A cell of atoms as a part of the whole system.
class CellList {
public:
    // Allocate memory for the cell list and return an empty configuration.
    CellList(const uint64_t capacity, const RVec origin, const RVec size);

    /************
    * Functions *
    *************/
    // Return the number of atoms in the list.
    uint64_t num_atoms(void) const { return natoms; };

    // Add an atom with input position to the system.
    void add_atom(const real x, const real y, const real z);

    Atom get_atom(const size_t index) const;

    /************
    * Variables *
    *************/
    std::vector<real> xs; // Positions (x), 1 elem per dimension: relative to the origin
    std::vector<real> vs; // Velocities (v)
    std::vector<real> fs; // Forces (f)
    std::vector<real> fs_prev; // Forces at previous step for Velocity Verlet

    // CellList origin.
    RVec origin;

    // CellList size.
    RVec size;

    // Neighbouring cell list indices in a collection (ie. a `System` struct)
    // which this list will interact with and add forces *to*. Since the force
    // calculation is symmetric we only have to calculate it once, going
    // *from* a list *to* another. Thus all list in this collection should
    // not have this list's index in its corresponding `to_neighbours` object.
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
    // The system is split into (maybe) several cell lists containing the atoms.
    //
    // If there are several cell lists, they should be non-overlapping
    // and added to this collection in the order of ZYX, meaning that
    // the indexing increases fastest along Z and slowest along X.
    //
    // The vector index of the cell list at (ix, iy, iz) thus is calculated as
    // i = ix * ny * nz + iy * nz + iz.
    std::vector<CellList> cell_lists;

    // System box size.
    RVec box_size;

    // System shape.
    IVec shape;

    // Size of cells in cell list.
    RVec cell_size;

    // Title of system.
    std::string title;
};

// Read a configuration from a Gromos formatted file.
System read_conf_from_grofile(const std::string& filename);

// Write the configuration to a Gromos formatted file.
void write_conf_to_grofile(const System& system, const std::string& path);

void create_cell_lists(System& system, const real rcut);

#endif // SYSTEM_CONF_H
