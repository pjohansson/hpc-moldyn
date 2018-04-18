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

// Add two vectors: r = r1 + r2.
RVec rvec_add(const RVec r1, const RVec r2);

// Subtract a vector from another: r = r1 - r2.
RVec rvec_sub(const RVec r1, const RVec r2);

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

    void update_num_atoms(void) {
        natoms = static_cast<uint64_t>(xs.size() / NDIM);
    }

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

// Read a configuration from a Gromos formatted file. Gromos is written
// with positions in nm, so use the sigma factor to convert.
System read_conf_from_grofile(const std::string& filename, const real sigma);

// Write the configuration to a Gromos formatted file. Gromos is written
// with positions in nm, so use the sigma factor to convert.
void write_conf_to_grofile(const System& system,
                           const std::string& path,
                           const real sigma);

// Construct the (3D) cell list configuration for the input system by splitting
// it into cells based on the input `target_size`. The cells will have sides
// that are *at least* that big in each dimension. Atoms are moved into their
// correct cells.
void create_cell_lists(System& system, const real target_size);

// Update the cell lists by moving atoms to their correct cells. For every
// atom, the position and velocity is moved as-is, the forces are set as
// the previous force of the atom, and the new current force is set to 0.
// This makes the update perfectly compatible with the Velocity Verlet scheme.
void update_cell_lists(System& system);

// Generate velocities for all atoms in the system around an input reference
// temperature.
void gen_system_velocities(System& system, const real Tref);

#endif // SYSTEM_CONF_H
