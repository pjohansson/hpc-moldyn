#include <mpi.h>
#include <ostream>
#include <chrono>
#include <thread>

#include "params.h"

#ifndef MPI_IMPL_H
#define MPI_IMPL_H

#ifdef USE_DOUBLE
#define MPI_MY_REAL_SIZE MPI_DOUBLE
#else
#define MPI_MY_REAL_SIZE MPI_FLOAT
#endif

constexpr size_t MASTER = 0;

/******************
 * DEBUG PRINTING *
 ******************/

#define MSDELAY 20

// Master prints
#define MASTER_PRINT(MPI_COMM, BODY) { \
    if (is_master(MPI_COMM)) \
    { \
        BODY \
    } \
    MPI_Barrier(MPI_COMM.comm_world); \
}

#define MASTER_PRINT_OBJECT(MPI_COMM, OBJECT) { \
    MASTER_PRINT(MPI_COMM, std::cerr << OBJECT << '\n';) \
}

// For every rank: Print the current rank number and evaluate the body.
#define MPI_RANK_PRINT(MPI_COMM, BODY) { \
    for (size_t RANK_ = 0; RANK_ < MPI_COMM.num_ranks; ++RANK_) \
    { \
        if (is_rank(RANK_, MPI_COMM)) \
        { \
            std::cerr << "rank " << RANK_ << ": "; \
            BODY \
            std::cerr << '\n' << std::flush; \
        } \
        std::this_thread::sleep_for(std::chrono::milliseconds(MSDELAY)); \
        MPI_Barrier(MPI_COMM.comm_world); \
    } \
}

// For every rank: Print the current rank number and a vector on it.
#define MPI_VEC_PRINT(MPI_COMM, VEC) { \
    MPI_RANK_PRINT(MPI_COMM, \
        std::cerr << '['; \
        for (const auto& VALUE_ : VEC) \
        { \
            std::cerr << ' ' << VALUE_; \
        } \
        std::cerr << " ]" << std::flush; \
    ) \
}

// For every rank: Print the current rank number and a vector of vectors on it.
#define MPI_VEC_VEC_PRINT(MPI_COMM, VEC) { \
    MPI_RANK_PRINT(MPI_COMM, \
        std::cerr << "\n[\n"; \
        for (const auto& INNER_VEC_ : VEC) \
        { \
            std::cerr << "  ["; \
            for (const auto& VALUE_ : INNER_VEC_) \
            { \
                std::cerr << ' ' << VALUE_;  \
            } \
            std::cerr << " ]\n"; \
        } \
        std::cerr << "]\n" << std::flush; \
    ) \
}

// For every rank: Print the current rank number and evaluate the body,
// without an MPI barrier.
#define MPI_RANK_PRINT_NOBARRIER(MPI_COMM, BODY) { \
    for (size_t RANK_ = 0; RANK_ < MPI_COMM.num_ranks; ++RANK_) \
    { \
        if (is_rank(RANK_, MPI_COMM)) \
        { \
            std::cerr << "rank " << RANK_ << ": "; \
            BODY \
            std::cerr << '\n' << std::flush; \
        } \
    } \
}

// For every rank: Print the current rank number and a vector on it.
#define MPI_VEC_PRINT_NOBARRIER(MPI_COMM, VEC) { \
    MPI_RANK_PRINT_NOBARRIER(MPI_COMM, \
        std::cerr << '['; \
        for (const auto& VALUE_ : VEC) \
        { \
            std::cerr << ' ' << VALUE_; \
        } \
        std::cerr << " ]" << std::flush; \
    ) \
}

struct MPICellComm {
    // Handle to the root of the communicator, which is that of
    // the rank that owns the cell.
    int root;

    // Handle to the communicator for the cell between its owning rank
    // and all those it will be shared with.
    MPI_Comm comm;

    // A list of the ranks the cell will be sent to.
    std::vector<size_t> to_ranks;
};

struct MPIRank {
    size_t rank,
           num_ranks;

    // For each cell in the system, the rank which owns it.
    std::vector<size_t> cell_parent_mpi_ranks;

    // For each MPI rank, which cells it owns. Every rank has this information
    // for every other rank, too.
    std::vector<std::vector<size_t>> mpi_rank_owned_cells;

    // For each MPI rank, which cells it does not own. The complement
    // to the owned field. Every rank has this information for every other
    // rank, too.
    std::vector<std::vector<size_t>> mpi_rank_non_owned_cells;

    // For each cell in the system: a communication record with all ranks
    // it will be shared with for the force calculation and a handle
    // to the owning rank inside the communication record.
    //
    // NOTE: The MPI handles are allocated on the heap *only for those processing
    // elements which are in the groups.* Be careful if freeing (right
    // now they'll all be freed at the end of runtime, since we don't
    // update them).
    std::vector<MPICellComm> cell_comm_groups;

    // For each MPI rank, which cells it will receive positions of for 
    // interaction calculations.
    std::vector<std::vector<size_t>> mpi_rank_received_cells;

    // For each MPI rank, which cells it will transmit positions of 
    // for interaction calculations. 
    std::vector<std::vector<size_t>> mpi_rank_sending_cells;

    // Communication record for all ranks participating in the simulation.
    // Kept separately from MPI_COMM_WORLD since some tests will not want 
    // to run on all ranks which necessitates specifying a subset.
    MPI_Comm comm_world;
};

bool is_master(const MPIRank& mpi_comm);
bool is_rank(const size_t rank, const MPIRank& mpi_comm);
bool init_MPI(MPIRank& mpi_comm);
bool sync_options(Options& opts, const MPIRank& mpi_comm);
bool sync_forcefield(ForceField& ff, const MPIRank& mpi_comm);

// Initial setup of the MPI system.
void mpi_init_cell_lists_and_transfer_atoms(System& system, MPIRank& mpi_comm);
void mpi_fill_communication_data(MPIRank& mpi_comm, const System& system);

// Interaction cell list synchronization
void mpi_synchronize_interaction_cell_lists(System& system,
                                            const MPIRank& mpi_comm);
void mpi_collect_forces_from_interaction_cell_lists(System& system,
                                                    const MPIRank& mpi_comm);
void reset_received_cell_lists(System& system, const MPIRank& mpi_comm);

// Move atoms to other ranks
void mpi_move_atoms_to_owning_ranks(System& system, const MPIRank& mpi_comm);
void mpi_collect_atoms_to_master(System& system, const MPIRank& mpi_comm);

#endif // MPI_IMPL_H
