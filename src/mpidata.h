#include <mpi.h>

#ifndef MPI_DATA_H
#define MPI_DATA_H

constexpr size_t MASTER = 0;

struct MPIRank {
    size_t rank,
           num_ranks;

    // set or vector of:
    // * owned cells
    // * information about which cells the rank needs for the integrator
    // * information about which other ranks own those cells
};

// function that takes the full system cell list from the master
// and divides the cells onto the ranks evenly

// etc

#endif // MPI_DATA_H
