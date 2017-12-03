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

constexpr char ATOM_NAME[2] = "C";
constexpr char RESIDUE_NAME[4] = "SOL";

class SystemConf {
public:
    // Allocate memory for the system and return an empty configuration.
    SystemConf(uint64_t capacity);

    /************
    * Functions *
    *************/
    // Return the number of atoms in the system.
    uint64_t num_atoms(void) const { return natoms; };

    // Add an atom with input position to the system.
    void add_atom(const real x, const real y, const real z);

    // Set the box size.
    void set_box(const real x, const real y, const real z);

    /************
    * Variables *
    *************/
    std::vector<real> xs; // Positions (x), 1 elem per dimension
    std::vector<real> vs; // Velocities (v)
    std::vector<real> fs; // Forces (f)

    // Box size.
    RVec box;

    // Title of system.
    std::string title;

private:
    uint64_t natoms;
};

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

    // Box origin.
    RVec origin;

    // Box size.
    RVec size;

private:
    uint64_t natoms;
};

// Read a configuration from a Gromos formatted file.
SystemConf read_conf_from_grofile(const std::string& filename);

// Write the configuration to a Gromos formatted file.
void write_conf_to_grofile(const SystemConf& conf, const std::string& path);

#endif // SYSTEM_CONF_H
