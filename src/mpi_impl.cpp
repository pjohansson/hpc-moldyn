#include <algorithm>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <set>
#include <vector>

#include "mpi_impl.h"

/*************
 * UTILITIES *
 *************/

bool is_master(const MPIRank& mpi_comm)
{
    return mpi_comm.rank == MASTER;
}

bool is_rank(const size_t rank, const MPIRank& mpi_comm)
{
    return mpi_comm.rank == rank;
}


/**************
 * BASE SETUP *
 **************/

bool init_MPI(MPIRank& mpi_comm)
{
    int rank, num_ranks;

    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    {
        return false;
    }

    if (MPI_Comm_size(MPI_COMM_WORLD, &num_ranks) != MPI_SUCCESS)
    {
        return false;
    }

    mpi_comm.rank = static_cast<size_t>(rank);
    mpi_comm.num_ranks = static_cast<size_t>(num_ranks);

    return true;
}

bool sync_options(Options& opts, const MPIRank& mpi_comm)
{
    MPI_Bcast(&opts.dt, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&opts.dt2, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&opts.gen_temp, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    MPI_Bcast(&opts.energy_calc, 1, MPI_UINT64_T, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&opts.num_steps, 1, MPI_UINT64_T, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&opts.traj_stride, 1, MPI_UINT64_T, MASTER, MPI_COMM_WORLD);

    auto gen_velocities = static_cast<uint64_t>(opts.gen_velocities);
    MPI_Bcast(&gen_velocities, 1, MPI_UINT64_T, MASTER, MPI_COMM_WORLD);
    opts.gen_velocities = gen_velocities;

    return true;
}


/**********************************************
 * STATIC BOOK KEEPING FOR MPI RANK OWNERSHIP *
 **********************************************/

// Get all the cell indices belonging to the input rank.
static std::vector<size_t> mpi_get_cell_indices(const uint64_t num_cells,
                                                const uint64_t rank,
                                                const uint64_t num_mpi_ranks)
{
    const auto max_rank_cells = static_cast<uint64_t>(ceil(
        static_cast<double>(num_cells) / static_cast<double>(num_mpi_ranks)
    ));

    auto num_ranks_with_max = static_cast<uint64_t>(num_cells % num_mpi_ranks);

    if (num_cells % max_rank_cells == 0)
    {
        num_ranks_with_max = num_mpi_ranks;
    }

    uint64_t first = 0,
             last = 0;

    if (rank < num_ranks_with_max)
    {
        first = rank * max_rank_cells;
        last = first + max_rank_cells;
    }
    else
    {
        first = num_ranks_with_max * max_rank_cells
            + (rank - num_ranks_with_max) * (max_rank_cells - 1);
        last = first + (max_rank_cells - 1);
    }

    std::vector<size_t> indices;

    for (auto index = first; index < last; ++index)
    {
        indices.push_back(static_cast<size_t>(index));
    }

    return indices;
}

// Given a list of owned cell for every MPI rank, return a list of rank
// indices for every cell.
static std::vector<size_t>
get_cell_mpi_ranks(const std::vector<std::vector<uint64_t>> rank_indices)
{
    std::vector<size_t> cell_ranks;

    size_t cell = 0;
    bool found = true;

    // Loop and look for the next cell index, until it is not found on any rank
    while (found)
    {
        found = false;
        uint64_t rank = 0;

        for (const auto indices : rank_indices)
        {
            const auto pos = std::find(indices.cbegin(), indices.cend(), cell);

            if (pos != indices.cend())
            {
                cell_ranks.push_back(rank);
                ++cell;
                found = true;

                break; // the for loop
            }

            ++rank;
        }
    }

    return cell_ranks;
}

// Return a list of the cells each rank owns.
static std::vector<std::vector<size_t>>
get_owned_cells_per_rank(const size_t num_cells, const size_t num_ranks)
{
    std::vector<std::vector<size_t>> owned_cells_per_rank;

    for (size_t rank = 0; rank < num_ranks; ++rank)
    {
        owned_cells_per_rank.push_back(mpi_get_cell_indices(
            num_cells, rank, num_ranks
        ));
    }

    return owned_cells_per_rank;
}

// Return a list of the cells which every rank does *not* own.
static std::vector<std::vector<size_t>> get_non_owned_cells_per_rank(
    const std::vector<std::vector<size_t>>& owned_cells_per_rank,
    const size_t                            num_cells
)
{
    std::vector<std::vector<size_t>> non_owned_cells_per_rank;

    std::set<size_t> all_cells;
    for (size_t i = 0; i < num_cells; ++i)
    {
        all_cells.insert(i);
    }

    for (const auto owned_cells : owned_cells_per_rank)
    {
        const std::set<size_t> owned_cells_set(
            owned_cells.cbegin(), owned_cells.cend()
        );

        std::set<size_t> complement;

        set_difference(
            all_cells.cbegin(), all_cells.cend(),
            owned_cells_set.cbegin(), owned_cells_set.cend(),
            std::inserter(complement, complement.begin())
        );

        std::vector<size_t> non_owned(complement.cbegin(), complement.cend());
        non_owned_cells_per_rank.push_back(non_owned);
    }

    return non_owned_cells_per_rank;
}

// Divide the system cells onto the present MPI ranks.
// Fill in the cells which every rank owns (`mpi_rank_owned_cells`)
// and which rank is the parent of each cell (`cell_parent_mpi_ranks`).
static void fill_mpi_rank_and_cell_ownership(MPIRank&     mpi_comm,
                                             const size_t num_cells)
{
    mpi_comm.mpi_rank_owned_cells = get_owned_cells_per_rank(
        num_cells, mpi_comm.num_ranks
    );
    mpi_comm.mpi_rank_non_owned_cells = get_non_owned_cells_per_rank(
        mpi_comm.mpi_rank_owned_cells, num_cells
    );
    mpi_comm.cell_parent_mpi_ranks = get_cell_mpi_ranks(
        mpi_comm.mpi_rank_owned_cells
    );
}

// For every cell, create an MPI_Comm group which it will be syncing from
// its parent, to every other rank that needs it.
static void mpi_create_cell_comm_groups(MPIRank& mpi_comm,
                                        const std::vector<CellList>& cell_lists)
{
    MPI_Group buf_group, world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    size_t i = 0;

    for (const auto& list : cell_lists)
    {
        // Initialize the collection with the parent of the cell,
        // since it will always be present. Then fill in the rest.
        const auto current_rank = static_cast<int>(
            mpi_comm.cell_parent_mpi_ranks.at(i)
        );
        std::set<size_t> to_ranks_set;

        for (const auto to_cell : list.to_neighbours)
        {
            const auto to_rank
                = mpi_comm.cell_parent_mpi_ranks.at(to_cell);

            if (static_cast<int>(to_rank) != current_rank)
            {
                to_ranks_set.insert(static_cast<size_t>(to_rank));
            }
        }

        const std::vector<size_t> to_ranks(
            to_ranks_set.cbegin(), to_ranks_set.cend()
        );

        std::vector<int> comm_group_ranks { static_cast<int>(current_rank) };
        for (const auto rank : to_ranks_set)
        {
            comm_group_ranks.push_back(static_cast<int>(rank));
        }

        MPI_Group_incl(
            world_group, comm_group_ranks.size(), comm_group_ranks.data(), &buf_group
        );

        MPI_Comm buf_comm;
        MPI_Comm_create(MPI_COMM_WORLD, buf_group, &buf_comm);

        // Get the cell root (owner) in the constructed communicator group
        // by translating from MPI_COMM_WORLD
        int cell_comm_root = -1;
        MPI_Group_translate_ranks(world_group, 1, &current_rank,
            buf_group, &cell_comm_root);

        MPICellComm cell_comm { cell_comm_root, buf_comm, to_ranks };
        mpi_comm.cell_comm_groups.push_back(cell_comm);

        MPI_Group_free(&buf_group);
        ++i;
    }

    MPI_Group_free(&world_group);
}

// Add information about which cells every rank will be sending and receiving
// to calculate the interactions.
static void mpi_create_sending_and_receiving_cell_lists_for_ranks(
    MPIRank& mpi_comm,
    const std::vector<CellList>& cell_lists
)
{
    // Using the to_neighbours list and list of parent ranks for every cell,
    // collect the cells which will be sent from each rank, and those which
    // will be received from other ranks, as sets to avoid double counting.
    // The indexing here is per MPI rank.
    std::vector<std::set<size_t>> received_cells_per_rank_set(
        mpi_comm.num_ranks, std::set<size_t> {}
    );
    std::vector<std::set<size_t>> sending_cells_per_rank_set(
        mpi_comm.num_ranks, std::set<size_t> {}
    );

    for (size_t i = 0; i < cell_lists.size(); ++i)
    {
        const auto& list = cell_lists.at(i);
        const auto parent_rank = mpi_comm.cell_parent_mpi_ranks.at(i);

        for (const auto to_cell : list.to_neighbours)
        {
            const auto to_rank = mpi_comm.cell_parent_mpi_ranks.at(to_cell);

            if (to_rank != parent_rank)
            {
                received_cells_per_rank_set.at(to_rank).insert(i);
                sending_cells_per_rank_set.at(parent_rank).insert(i);
            }
        }
    }

    std::vector<std::vector<size_t>> received_cells_per_rank;
    std::vector<std::vector<size_t>> sending_cells_per_rank;

    for (const auto& received_set : received_cells_per_rank_set)
    {
        const std::vector<size_t> received_cells(
            received_set.cbegin(), received_set.cend()
        );
        received_cells_per_rank.push_back(received_cells);
    }

    for (const auto& sending_set : sending_cells_per_rank_set)
    {
        const std::vector<size_t> sending_cells(
            sending_set.cbegin(), sending_set.cend()
        );
        sending_cells_per_rank.push_back(sending_cells);
    }

    mpi_comm.mpi_rank_received_cells = received_cells_per_rank;
    mpi_comm.mpi_rank_sending_cells = sending_cells_per_rank;
}


/**************************
 * SYSTEM CELL LIST SETUP *
 **************************/

// Given a system in which rank 0 has all the cell lists, create
// same structure on every other rank.
static void mpi_create_all_cell_lists(System& system,
                                      const uint64_t num_cells,
                                      const MPIRank& mpi_comm)
{
    RVec cell_size,
         origin;

    if (is_master(mpi_comm))
    {
        cell_size = system.cell_lists.at(0).size;
    }

    MPI_Bcast(cell_size.data(), NDIM, MPI_MY_REAL_SIZE, MASTER, MPI_COMM_WORLD);

    for (size_t index = 0; index < num_cells; ++index)
    {
        if (is_master(mpi_comm))
        {
            origin = system.cell_lists.at(index).origin;
        }

        MPI_Bcast(
            origin.data(), NDIM, MPI_MY_REAL_SIZE, MASTER, MPI_COMM_WORLD
        );

        if (!is_master(mpi_comm))
        {
            CellList list {0, origin, cell_size};
            system.cell_lists.push_back(list);
        }
    }
}

// Given a system in which rank 0 owns all the cell lists, transfer them
// to their parent ranks.
void mpi_init_cell_lists_and_transfer(System& system, MPIRank& mpi_comm)
{
    auto num_cells = static_cast<uint64_t>(system.cell_lists.size());
    MPI_Bcast(&num_cells, 1, MPI_UINT64_T, MASTER, MPI_COMM_WORLD);

    fill_mpi_rank_and_cell_ownership(mpi_comm, num_cells);
    mpi_create_all_cell_lists(system, num_cells, mpi_comm);
    mpi_move_atoms_to_owning_ranks(system, mpi_comm);

    if (is_master(mpi_comm))
    {
        const auto& non_owned = mpi_comm.mpi_rank_non_owned_cells.at(MASTER);
        for (const auto& index_cell : non_owned)
        {
            system.cell_lists.at(index_cell).resize_atom_list(0);
        }
    }
}


/***************************************************
 * MPI TRANSFER OF INTERACTION CELLS TO NEIGHBOURS *
 ***************************************************/

static std::vector<MPI_Request> mpi_send_numbers(
    const std::vector<CellList>& cell_lists,
    const MPIRank&               mpi_comm
)
{
    const auto& send_cells = mpi_comm.mpi_rank_sending_cells.at(mpi_comm.rank);
    std::vector<MPI_Request> mpi_send_requests(send_cells.size(), nullptr);

    for (size_t i = 0; i < send_cells.size(); ++i)
    {
        const auto& index_cell = send_cells.at(i);
        auto& request = mpi_send_requests.at(i);

        const auto num = static_cast<uint64_t>(
            cell_lists.at(index_cell).num_atoms()
        );
        const auto to_ranks = mpi_comm.cell_comm_groups.at(index_cell).to_ranks;

        for (const auto rank : to_ranks)
        {
            MPI_Isend(
                &num, 1, MPI_UINT64_T,
                rank, index_cell, MPI_COMM_WORLD, &request
            );
        }
    }

    return mpi_send_requests;
}

static std::vector<MPI_Request> mpi_recv_numbers(
    std::vector<size_t>&         num_recv_per_cell,
    const std::vector<CellList>& cell_lists,
    const MPIRank&               mpi_comm
)
{
    const auto& recv_cells = mpi_comm.mpi_rank_received_cells.at(mpi_comm.rank);
    std::vector<MPI_Request> mpi_recv_requests(recv_cells.size(), nullptr);

    for (size_t i = 0; i < recv_cells.size(); ++i)
    {
        const auto& index_cell = recv_cells.at(i);
        auto& request = mpi_recv_requests.at(i);

        const auto from_rank = mpi_comm.cell_parent_mpi_ranks.at(index_cell);

        MPI_Irecv(&num_recv_per_cell.at(i), 1, MPI_UINT64_T,
            from_rank, index_cell, MPI_COMM_WORLD, &request);
    }

    return mpi_recv_requests;
}

// For every individual MPI rank, return the number of positions
// that will be transferred from other ranks for every cell in
// the receiving cells list.
static std::vector<size_t> mpi_sync_number_of_transmitted_atoms(
    const std::vector<CellList>& cell_lists,
    const MPIRank&               mpi_comm
)
{
    const auto& recv_cells = mpi_comm.mpi_rank_received_cells.at(mpi_comm.rank);
    std::vector<size_t> num_recv_per_cell(recv_cells.size(), 0);

    auto mpi_send_requests = mpi_send_numbers(cell_lists, mpi_comm);
    auto mpi_recv_requests = mpi_recv_numbers(
        num_recv_per_cell, cell_lists, mpi_comm
    );

    // MPI_Waitall will free the requests, so we don't have to.
    MPI_Waitall(
        static_cast<int>(mpi_send_requests.size()),
        mpi_send_requests.data(),
        MPI_STATUSES_IGNORE
    );
    MPI_Waitall(
        static_cast<int>(mpi_recv_requests.size()),
        mpi_recv_requests.data(),
        MPI_STATUSES_IGNORE
    );

    MPI_Barrier(MPI_COMM_WORLD);

    return num_recv_per_cell;
}

static void reserve_memory_for_received_lists(
    std::vector<CellList>&    cell_lists,
    const std::vector<size_t> num_recv_per_cell,
    const std::vector<size_t> recv_cell_indices
)
{
    for (size_t i = 0; i < num_recv_per_cell.size(); ++i)
    {
        const auto& num = num_recv_per_cell.at(i);
        const auto& index_cell = recv_cell_indices.at(i);

        cell_lists.at(index_cell).resize_atom_list_force_calc(num);
    }
}

static std::vector<MPI_Request> mpi_send_cell_list_positions(
    std::vector<CellList>& cell_lists,
    const MPIRank&         mpi_comm
)
{
    const auto& send_cells = mpi_comm.mpi_rank_sending_cells.at(mpi_comm.rank);
    std::vector<MPI_Request> mpi_requests(send_cells.size(), nullptr);

    for (size_t i = 0; i < send_cells.size(); ++i)
    {
        const auto& index_cell = send_cells.at(i);
        auto& request = mpi_requests.at(i);

        auto& list = cell_lists.at(index_cell);
        const auto& comm_group = mpi_comm.cell_comm_groups.at(index_cell).comm;
        const auto& comm_root = mpi_comm.cell_comm_groups.at(index_cell).root;

        MPI_Ibcast(
            list.xs.data(), list.xs.size(), MPI_MY_REAL_SIZE,
            comm_root, comm_group, &request
        );
    }

    return mpi_requests;
}

static std::vector<MPI_Request> mpi_recv_cell_list_positions(
    std::vector<CellList>& cell_lists,
    const MPIRank&         mpi_comm
)
{
    const auto& recv_cells = mpi_comm.mpi_rank_received_cells.at(mpi_comm.rank);
    std::vector<MPI_Request> mpi_requests(recv_cells.size(), nullptr);

    for (size_t i = 0; i < recv_cells.size(); ++i)
    {
        const auto& index_cell = recv_cells.at(i);
        auto& request = mpi_requests.at(i);

        auto& list = cell_lists.at(index_cell);
        const auto& comm_group = mpi_comm.cell_comm_groups.at(index_cell).comm;
        const auto& comm_root = mpi_comm.cell_comm_groups.at(index_cell).root;

        MPI_Ibcast(
            list.xs.data(), list.xs.size(), MPI_MY_REAL_SIZE,
            comm_root, comm_group, &request
        );
    }

    return mpi_requests;
}

static void mpi_transmit_interaction_cell_lists(
    std::vector<CellList>& cell_lists,
    const MPIRank&         mpi_comm
)
{
    auto mpi_send_requests = mpi_send_cell_list_positions(cell_lists, mpi_comm);
    auto mpi_recv_requests = mpi_recv_cell_list_positions(cell_lists, mpi_comm);

    MPI_Waitall(
        static_cast<int>(mpi_send_requests.size()),
        mpi_send_requests.data(),
        MPI_STATUSES_IGNORE
    );
    MPI_Waitall(
        static_cast<int>(mpi_recv_requests.size()),
        mpi_recv_requests.data(),
        MPI_STATUSES_IGNORE
    );
    MPI_Barrier(MPI_COMM_WORLD);
}

static std::vector<MPI_Request> mpi_send_cell_list_forces(
    std::vector<CellList>& cell_lists,
    const MPIRank&         mpi_comm
)
{
    // Transmit the forces *from the received cells*, not the sending
    const auto& send_cells = mpi_comm.mpi_rank_received_cells.at(mpi_comm.rank);
    std::vector<MPI_Request> mpi_requests(send_cells.size(), nullptr);

    for (size_t i = 0; i < send_cells.size(); ++i)
    {
        const auto& index_cell = send_cells.at(i);
        auto& request = mpi_requests.at(i);

        auto& list = cell_lists.at(index_cell);
        const auto& comm_group = mpi_comm.cell_comm_groups.at(index_cell).comm;
        const auto& comm_root = mpi_comm.cell_comm_groups.at(index_cell).root;

        MPI_Ireduce(
            list.fs.data(), nullptr, list.fs.size(), MPI_MY_REAL_SIZE,
            MPI_SUM, comm_root, comm_group, &request
        );
    }

    return mpi_requests;
}

static std::vector<MPI_Request> mpi_recv_cell_list_forces(
    std::vector<CellList>& cell_lists,
    const MPIRank&         mpi_comm
)
{
    // Receive the forces *in the sending cells*, not the receiving
    const auto& recv_cells = mpi_comm.mpi_rank_sending_cells.at(mpi_comm.rank);
    std::vector<MPI_Request> mpi_requests(recv_cells.size(), nullptr);

    for (size_t i = 0; i < recv_cells.size(); ++i)
    {
        const auto& index_cell = recv_cells.at(i);
        auto& request = mpi_requests.at(i);

        auto& list = cell_lists.at(index_cell);
        const auto& comm_group = mpi_comm.cell_comm_groups.at(index_cell).comm;
        const auto& comm_root = mpi_comm.cell_comm_groups.at(index_cell).root;

        // Use MPI_IN_PLACE to use the send-buffer on the root rank
        // (which will be these) as the recv-buffer.
        MPI_Ireduce(
            MPI_IN_PLACE, list.fs.data(), list.fs.size(), MPI_MY_REAL_SIZE,
            MPI_SUM, comm_root, comm_group, &request
        );
    }

    return mpi_requests;
}

void mpi_synchronize_interaction_cell_lists(System& system,
                                            const MPIRank& mpi_comm)
{
    auto& cell_lists = system.cell_lists;

    const auto num_recv_per_cell = mpi_sync_number_of_transmitted_atoms(
        cell_lists, mpi_comm
    );

    const auto& rank_received_cells
        = mpi_comm.mpi_rank_received_cells.at(mpi_comm.rank);

    // DEBUG: Assert that things are going well (?)
    if (num_recv_per_cell.size() != rank_received_cells.size())
    {
        MPI_RANK_PRINT(mpi_comm,
            std::cerr
                << "WARNING: when transmitting cell lists for "
                   "the interactions, the rank is trying to receive "
                   "positions in more cells (" << num_recv_per_cell.size()
                << ") than it was expecting (" << rank_received_cells.size()
                << "). Something is wrong with the MPI book keeping.\n";
        )

        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    reserve_memory_for_received_lists(
        cell_lists, num_recv_per_cell, rank_received_cells
    );

    mpi_transmit_interaction_cell_lists(cell_lists, mpi_comm);
}

void mpi_collect_forces_from_interaction_cell_lists(System& system,
                                                    const MPIRank& mpi_comm)
{
    auto& cell_lists = system.cell_lists;

    auto mpi_send_requests = mpi_send_cell_list_forces(cell_lists, mpi_comm);

    MPI_Waitall(
        static_cast<int>(mpi_send_requests.size()),
        mpi_send_requests.data(),
        MPI_STATUSES_IGNORE
    );

    auto mpi_recv_requests = mpi_recv_cell_list_forces(cell_lists, mpi_comm);

    MPI_Waitall(
        static_cast<int>(mpi_recv_requests.size()),
        mpi_recv_requests.data(),
        MPI_STATUSES_IGNORE
    );

    MPI_Barrier(MPI_COMM_WORLD);
}

// For the calling MPI rank:
// Clear all the cells that receive atoms from other ranks.
void reset_received_cell_lists(System& system, const MPIRank& mpi_comm)
{
    auto& cell_lists = system.cell_lists;

    const auto received_cell_indices
        = mpi_comm.mpi_rank_received_cells.at(mpi_comm.rank);

    for (const auto index_cell : received_cell_indices)
    {
        cell_lists.at(index_cell).resize_atom_list(0);
    }
}


/**************************************
 * MPI TRANSFER OF ATOMS TO NEW RANKS *
 **************************************/

// For the calling MPI rank:
// For every rank, get the cells to transmit to them.
static std::vector<std::vector<uint64_t>>
get_cells_to_transmit(const std::vector<CellList>& cell_lists,
                      const MPIRank&               mpi_comm)
{
    const auto& non_owned_cells_per_rank = mpi_comm.mpi_rank_non_owned_cells.at(
        mpi_comm.rank
    );
    std::vector<std::vector<uint64_t>> cells_to_transmit_per_rank(
        mpi_comm.num_ranks
    );

    for (const auto index_cell : non_owned_cells_per_rank)
    {
        const auto& list = cell_lists.at(index_cell);

        if (list.num_atoms() > 0)
        {
            const auto& to_rank = mpi_comm.cell_parent_mpi_ranks.at(index_cell);
            cells_to_transmit_per_rank.at(to_rank).push_back(index_cell);
        }
    }

    return cells_to_transmit_per_rank;
}

// Synchronize the number of cells that will be sent from each MPI rank
// to the current caller.
static std::vector<size_t> mpi_sync_number_of_transmitted_cells(
    const std::vector<std::vector<uint64_t>>& cells_to_transmit_per_rank,
    const MPIRank&                            mpi_comm
)
{
    std::vector<uint64_t> num_cells_to_transmit_per_rank(mpi_comm.num_ranks, 0);
    for (size_t rank = 0; rank < mpi_comm.num_ranks; ++rank)
    {
        num_cells_to_transmit_per_rank.at(rank) = static_cast<uint64_t>(
            cells_to_transmit_per_rank.at(rank).size()
        );
    }

    std::vector<uint64_t> num_cells_to_receive_per_rank(mpi_comm.num_ranks, 0);
    std::vector<MPI_Request> mpi_requests(mpi_comm.num_ranks, nullptr);

    for (size_t rank = 0; rank < mpi_comm.num_ranks; ++rank)
    {
        MPI_Iscatter(
            num_cells_to_transmit_per_rank.data(), 1, MPI_UINT64_T,
            num_cells_to_receive_per_rank.data() + rank, 1, MPI_UINT64_T,
            rank, MPI_COMM_WORLD, &mpi_requests.at(rank)
        );
    }

    MPI_Waitall(
        static_cast<int>(mpi_requests.size()),
        mpi_requests.data(),
        MPI_STATUSES_IGNORE
    );

    return num_cells_to_receive_per_rank;
}

struct CellRecvInfo {
    std::vector<uint64_t> cell_indices, num_atoms;
};

static MPI_Request
mpi_send_cell_indices_and_num_atoms(const std::vector<CellList>& cell_lists,
                                    const std::vector<uint64_t>& transmit_cells,
                                    const size_t                 to_rank)
{
    MPI_Request cells_request = nullptr,
                nums_request = nullptr;

    const auto num_cells = static_cast<int>(transmit_cells.size());

    // Collect the number of atoms every cell will transfer
    std::vector<uint64_t> transmit_nums(num_cells, 0);

    for (auto i = 0; i < num_cells; ++i)
    {
        const auto index_cell = transmit_cells.at(i);
        transmit_nums.at(i) = cell_lists.at(index_cell).num_atoms();
    }

    MPI_Isend(
        transmit_nums.data(), num_cells,
        MPI_UINT64_T, to_rank, 0,
        MPI_COMM_WORLD, &nums_request
    );

    MPI_Isend(
        transmit_cells.data(), num_cells,
        MPI_UINT64_T, to_rank, 1,
        MPI_COMM_WORLD, &cells_request
    );

    // Wait for the atom number send to be complete before returning,
    // since we don't want to deallocate the data before that.
    MPI_Wait(&nums_request, MPI_STATUS_IGNORE);

    return cells_request;
}

static std::vector<CellRecvInfo>
mpi_sync_moving_cell_information(
    const std::vector<CellList>&              cell_lists,
    const std::vector<std::vector<uint64_t>>& cells_to_transmit_per_rank,
    const std::vector<uint64_t>&              num_cells_to_receive_per_rank,
    const MPIRank&                            mpi_comm
)
{
    // Receive the cell indices and their number of atoms in this collection
    std::vector<CellRecvInfo> recv_cell_info_per_rank(mpi_comm.num_ranks);

    std::vector<MPI_Request> mpi_send_cells_requests(
        mpi_comm.num_ranks, MPI_REQUEST_NULL
    );
    std::vector<MPI_Request> mpi_recv_cells_requests(
        mpi_comm.num_ranks, MPI_REQUEST_NULL
    );
    std::vector<MPI_Request> mpi_recv_num_atoms_requests(
        mpi_comm.num_ranks, MPI_REQUEST_NULL
    );

    for (size_t rank = 0; rank < mpi_comm.num_ranks; ++rank)
    {
        if (!is_rank(rank, mpi_comm))
        {
            auto& info = recv_cell_info_per_rank.at(rank);
            const auto num = num_cells_to_receive_per_rank.at(rank);

            info.cell_indices.resize(num);
            info.num_atoms.resize(num);

            // Collect the cell indices which will be transmitted, and the
            // number of atoms for those cells
            const auto& transmit_cells = cells_to_transmit_per_rank.at(rank);

            if (transmit_cells.size() > 0)
            {
                mpi_send_cells_requests.at(rank)
                    = mpi_send_cell_indices_and_num_atoms(
                        cell_lists, transmit_cells, rank);
            }

            if (num > 0)
            {
                MPI_Irecv(
                    info.num_atoms.data(),
                    static_cast<int>(num), MPI_UINT64_T,
                    rank, 0,
                    MPI_COMM_WORLD, &mpi_recv_num_atoms_requests.at(rank)
                );

                MPI_Irecv(
                    info.cell_indices.data(),
                    static_cast<int>(num), MPI_UINT64_T,
                    rank, 1,
                    MPI_COMM_WORLD, &mpi_recv_cells_requests.at(rank)
                );
            }
        }
    }

    MPI_Waitall(
        static_cast<int>(mpi_comm.num_ranks),
        mpi_send_cells_requests.data(),
        MPI_STATUSES_IGNORE
    );
    MPI_Waitall(
        static_cast<int>(mpi_comm.num_ranks),
        mpi_recv_cells_requests.data(),
        MPI_STATUSES_IGNORE
    );
    MPI_Waitall(
        static_cast<int>(mpi_comm.num_ranks),
        mpi_recv_num_atoms_requests.data(),
        MPI_STATUSES_IGNORE
    );
    MPI_Barrier(MPI_COMM_WORLD);

    return recv_cell_info_per_rank;
}

// Reserve memory for the received atoms and return a list
// of the current number of atoms per cell for indexing purposes.
static std::vector<size_t> reserve_memory_for_received_atoms(
    std::vector<CellList>& cell_lists,
    const std::vector<CellRecvInfo>& recv_cell_info_per_rank
)
{
    std::vector<size_t> received_atoms_per_cell(cell_lists.size(), 0);

    for (const auto& info : recv_cell_info_per_rank)
    {
        for (size_t i = 0; i < info.num_atoms.size(); ++i)
        {
            received_atoms_per_cell.at(info.cell_indices.at(i))
                += info.num_atoms.at(i);
        }
    }

    std::vector<size_t> num_current_atoms_per_cell(cell_lists.size(), 0);

    for (size_t i = 0; i < cell_lists.size(); ++i)
    {
        auto& list = cell_lists.at(i);
        num_current_atoms_per_cell.at(i) = list.num_atoms();

        const auto new_num_atoms
            = list.num_atoms() + received_atoms_per_cell.at(i);

        list.resize_atom_list(new_num_atoms);
    }

    return num_current_atoms_per_cell;
}

struct SendRecvRequests {
    std::vector<MPI_Request> positions, velocities;
};

static SendRecvRequests mpi_send_atoms(
    std::vector<CellList>&                  cell_lists,
    const std::vector<std::vector<size_t>>& cells_to_transmit_per_rank,
    const MPIRank&                          mpi_comm
)
{
    std::vector<MPI_Request> mpi_pos_requests, mpi_vel_requests;

    for (size_t rank = 0; rank < mpi_comm.num_ranks; ++rank)
    {
        const auto& transmit_cells = cells_to_transmit_per_rank.at(rank);

        for (const auto& index_cell : transmit_cells)
        {
            auto& list = cell_lists.at(index_cell);

            MPI_Request pos_request = nullptr,
                        vel_request = nullptr;

            MPI_Isend(
                list.xs.data(), static_cast<int>(list.xs.size()), MPI_MY_REAL_SIZE,
                rank, index_cell,
                MPI_COMM_WORLD, &pos_request
            );

            MPI_Isend(
                list.vs.data(), static_cast<int>(list.vs.size()), MPI_MY_REAL_SIZE,
                rank, index_cell + cell_lists.size(),
                MPI_COMM_WORLD, &vel_request
            );

            mpi_pos_requests.push_back(pos_request);
            mpi_vel_requests.push_back(vel_request);
        }
    }

    SendRecvRequests send_requests {
        std::move(mpi_pos_requests),
        std::move(mpi_vel_requests)
    };

    return send_requests;
}

static SendRecvRequests mpi_receive_atoms(
    std::vector<CellList>&           cell_lists,
    const std::vector<CellRecvInfo>& recv_cell_info_per_rank,
    const MPIRank&                   mpi_comm
)
{
    auto num_atoms_per_cell = reserve_memory_for_received_atoms(
        cell_lists, recv_cell_info_per_rank
    );

    // For every rank, go through all their cells and receive both their
    // positions and velocities. Ensure that they are placed at the correct
    // place in every array by keeping track of the number received in
    // every cell.
    std::vector<MPI_Request> mpi_pos_requests, mpi_vel_requests;

    for (size_t rank = 0; rank < mpi_comm.num_ranks; ++rank)
    {
        const auto& info = recv_cell_info_per_rank.at(rank);
        const auto& num_cells_from_rank = info.cell_indices.size();

        for (size_t i = 0; i < num_cells_from_rank; ++i)
        {
            const auto& index_cell = info.cell_indices.at(i);
            const auto num_atoms = static_cast<int>(info.num_atoms.at(i));

            auto& list = cell_lists.at(index_cell);
            auto& num_atoms_present = num_atoms_per_cell.at(index_cell);
            const auto begin = NDIM * num_atoms_present;

            MPI_Request pos_request = nullptr,
                        vel_request = nullptr;

            // The message tag for the positions from each rank is
            // the cell index, the velocity tag is that plus the number
            // of cells in the system to ensure that it is unique.
            MPI_Irecv(
                list.xs.data() + begin, NDIM * num_atoms, MPI_MY_REAL_SIZE,
                rank, index_cell,
                MPI_COMM_WORLD, &pos_request
            );

            MPI_Irecv(
                list.vs.data() + begin, NDIM * num_atoms, MPI_MY_REAL_SIZE,
                rank, cell_lists.size() + index_cell,
                MPI_COMM_WORLD, &vel_request
            );

            num_atoms_present += static_cast<size_t>(num_atoms);

            mpi_pos_requests.push_back(pos_request);
            mpi_vel_requests.push_back(vel_request);
        }
    }

    SendRecvRequests recv_requests {
        std::move(mpi_pos_requests),
        std::move(mpi_vel_requests)
    };

    return recv_requests;
}

static void
mpi_move_atoms(
    std::vector<CellList>&                  cell_lists,
    const std::vector<std::vector<size_t>>& cells_to_transmit_per_rank,
    const std::vector<CellRecvInfo>&        recv_cell_info_per_rank,
    const MPIRank&                          mpi_comm
)
{
    auto mpi_send_requests = mpi_send_atoms(
        cell_lists, cells_to_transmit_per_rank, mpi_comm
    );

    auto mpi_recv_requests = mpi_receive_atoms(
        cell_lists, recv_cell_info_per_rank, mpi_comm
    );

    MPI_Waitall(
        static_cast<int>(mpi_send_requests.positions.size()),
        mpi_send_requests.positions.data(),
        MPI_STATUSES_IGNORE
    );
    MPI_Waitall(
        static_cast<int>(mpi_send_requests.velocities.size()),
        mpi_send_requests.velocities.data(),
        MPI_STATUSES_IGNORE
    );
    MPI_Waitall(
        static_cast<int>(mpi_recv_requests.positions.size()),
        mpi_recv_requests.positions.data(),
        MPI_STATUSES_IGNORE
    );
    MPI_Waitall(
        static_cast<int>(mpi_recv_requests.velocities.size()),
        mpi_recv_requests.velocities.data(),
        MPI_STATUSES_IGNORE
    );

    MPI_Barrier(MPI_COMM_WORLD);
}

void mpi_move_atoms_to_owning_ranks(System& system, const MPIRank& mpi_comm)
{
    const auto cells_to_transmit_per_rank = get_cells_to_transmit(
        system.cell_lists, mpi_comm
    );

    const auto num_cells_to_receive_per_rank
        = mpi_sync_number_of_transmitted_cells(
            cells_to_transmit_per_rank, mpi_comm);

    const auto recv_cell_info_per_rank = mpi_sync_moving_cell_information(
        system.cell_lists, cells_to_transmit_per_rank,
        num_cells_to_receive_per_rank, mpi_comm
    );

    mpi_move_atoms(
        system.cell_lists, cells_to_transmit_per_rank,
        recv_cell_info_per_rank, mpi_comm
    );
}

void mpi_collect_atoms_to_master(System& system, const MPIRank& mpi_comm)
{
    // Use the same framework as mpi_move_atoms_to_owning_ranks:
    // Get a list of cells to transmit to every rank, but only populate
    // the list with cells to send to master.
    std::vector<std::vector<uint64_t>> cells_to_transmit_per_rank(
        mpi_comm.num_ranks
    );

    if (!is_master(mpi_comm))
    {
        const auto& my_owned_cells
            = mpi_comm.mpi_rank_owned_cells.at(mpi_comm.rank);

        for (const auto& index_cell : my_owned_cells)
        {
            const auto& list = system.cell_lists.at(index_cell);

            if (list.num_atoms() > 0)
            {
                cells_to_transmit_per_rank.at(MASTER).push_back(index_cell);
            }
        }
    }

    const auto num_cells_to_receive_per_rank
        = mpi_sync_number_of_transmitted_cells(
            cells_to_transmit_per_rank, mpi_comm);

    const auto recv_cell_info_per_rank = mpi_sync_moving_cell_information(
        system.cell_lists, cells_to_transmit_per_rank,
        num_cells_to_receive_per_rank, mpi_comm
    );

    mpi_move_atoms(
        system.cell_lists, cells_to_transmit_per_rank,
        recv_cell_info_per_rank, mpi_comm
    );
}
