#include "tests/utils.h"

#include <iostream>

#include "src/mpi_impl.cpp"

using namespace std;

template<typename T>
static bool in_set(const T needle, const set<T>& haystack)
{
    return static_cast<bool>(haystack.find(needle) != haystack.cend());
}

template<typename T>
static bool in_set(const T needle, const vector<T>& haystack)
{
    return static_cast<bool>(
        find(haystack.cbegin(), haystack.cend(), needle) != haystack.cend()
    );
}

ADD_TEST(test_cell_indices_are_correctly_divided,
    const auto num_mpi_ranks = 4;
    const auto num_cells_A = 15;

    vector<uint64_t> rank0_indices { 0, 1, 2, 3 };
    ASSERT_EQ_VEC(mpi_get_cell_indices(num_cells_A, 0, num_mpi_ranks),
                  rank0_indices,
                  "rank 0 does not get correct indices for 15");

    vector<uint64_t> rank1_indices { 4, 5, 6, 7 };
    ASSERT_EQ_VEC(mpi_get_cell_indices(num_cells_A, 1, num_mpi_ranks),
                  rank1_indices,
                  "rank 1 does not get correct indices for 15");

    vector<uint64_t> rank2_indices { 8, 9, 10, 11 };
    ASSERT_EQ_VEC(mpi_get_cell_indices(num_cells_A, 2, num_mpi_ranks),
                  rank2_indices,
                  "rank 2 does not get correct indices for 15");

    vector<uint64_t> rank3_indices { 12, 13, 14 };
    ASSERT_EQ_VEC(mpi_get_cell_indices(num_cells_A, 3, num_mpi_ranks),
                  rank3_indices,
                  "rank 3 does not get correct indices for 15");

    const auto num_cells_B = 13;

    rank0_indices = vector<uint64_t> { 0, 1, 2, 3 };
    ASSERT_EQ_VEC(mpi_get_cell_indices(num_cells_B, 0, num_mpi_ranks),
                  rank0_indices,
                  "rank 0 does not get correct indices for 13");

    rank1_indices = vector<uint64_t> { 4, 5, 6 };
    ASSERT_EQ_VEC(mpi_get_cell_indices(num_cells_B, 1, num_mpi_ranks),
                  rank1_indices,
                  "rank 1 does not get correct indices for 13");

    rank2_indices = vector<uint64_t> { 7, 8, 9 };
    ASSERT_EQ_VEC(mpi_get_cell_indices(num_cells_B, 2, num_mpi_ranks),
                  rank2_indices,
                  "rank 2 does not get correct indices for 13");

    rank3_indices = vector<uint64_t> { 10, 11, 12 };
    ASSERT_EQ_VEC(mpi_get_cell_indices(num_cells_B, 3, num_mpi_ranks),
                  rank3_indices,
                  "rank 3 does not get correct indices for 13");
)

ADD_TEST(test_get_owned_cells_per_rank,
    const size_t num_ranks = 3;
    const size_t num_cells = 12;

    vector<vector<size_t>> expected_cells_per_rank;
    for (size_t rank = 0; rank < num_ranks; ++rank)
    {
        const auto indices = mpi_get_cell_indices(num_cells, rank, num_ranks);
        expected_cells_per_rank.push_back(indices);
    }

    const auto owned_cells_per_rank = get_owned_cells_per_rank(num_cells, num_ranks);

    for (size_t rank = 0; rank < num_ranks; ++rank)
    {
        ASSERT_EQ_VEC(
            owned_cells_per_rank.at(rank), expected_cells_per_rank.at(rank),
            "a ranks got the wrong owned cells"
        );
    }
)

ADD_TEST(test_get_owned_cells_per_rank_even_dividers,
    const size_t num_ranks = 4;
    const size_t num_cells_1 = 4;

    const auto owned_cells_per_rank_1 = get_owned_cells_per_rank(
        num_cells_1, num_ranks
    );

    for (size_t i = 0; i < 4; ++i)
    {
        ASSERT_EQ(
            owned_cells_per_rank_1.at(i).size(), 1,
            "a rank got the wrong number of cells"
        );
        ASSERT_EQ(owned_cells_per_rank_1.at(i)[0], i,
            "a rank got the wrong cell");
    }

    const size_t num_cells_2 = 8;

    const auto owned_cells_per_rank_2 = get_owned_cells_per_rank(
        num_cells_2, num_ranks
    );

    const vector<size_t> expected {0, 1, 2, 3, 4, 5, 6, 7};

    for (size_t i = 0; i < 4; ++i)
    {
        const vector<size_t> expected_slice (
            expected.cbegin() + 2 * i, expected.cbegin() + 2 * (i + 1)
        );

        ASSERT_EQ(
            owned_cells_per_rank_2.at(i).size(), 2,
            "a rank got the wrong number of cells"
        );
        ASSERT_EQ_VEC(owned_cells_per_rank_2.at(i), expected_slice,
            "a rank got the wrong cells");
    }
)

ADD_TEST(test_get_owned_cells_per_rank_odd_dividers,
    const size_t num_ranks = 3;
    const size_t num_cells = 320;

    const auto owned_cells_per_rank = get_owned_cells_per_rank(
        num_cells, num_ranks
    );

    vector<size_t> rank0_expected; // should get 107 cells: 0-106
    vector<size_t> rank1_expected; // 107 cells: 107-213
    vector<size_t> rank2_expected; // 106 cells: 214-319

    for (size_t i = 0; i <= 106; ++i)
    {
        rank0_expected.push_back(i);
    }
    for (size_t i = 107; i <= 213; ++i)
    {
        rank1_expected.push_back(i);
    }
    for (size_t i = 214; i <= 319; ++i)
    {
        rank2_expected.push_back(i);
    }

    ASSERT_EQ(owned_cells_per_rank.size(), 3,
        "incorrect number of ranks returned");
    ASSERT_EQ_VEC(owned_cells_per_rank.at(0), rank0_expected,
        "rank 0 does not have the correct indices");
    ASSERT_EQ_VEC(owned_cells_per_rank.at(1), rank1_expected,
        "rank 1 does not have the correct indices");
    ASSERT_EQ_VEC(owned_cells_per_rank.at(2), rank2_expected,
        "rank 2 does not have the correct indices");
)

ADD_TEST(test_get_cell_ranks_from_owned_cells_per_rank_lists,
    vector<vector<uint64_t>> owned_cells_per_rank;

    owned_cells_per_rank.push_back({0, 3, 6});
    owned_cells_per_rank.push_back({1, 4, 7});
    owned_cells_per_rank.push_back({2, 5});

    const auto cell_ranks = get_cell_mpi_ranks(owned_cells_per_rank);

    const vector<uint64_t> expected_parents { 0, 1, 2, 0, 1, 2, 0, 1 };
    ASSERT_EQ_VEC(cell_ranks, expected_parents, "cell ranks do not match the list");
)

ADD_TEST(test_get_non_owned_cells_for_ranks,
    const size_t num_cells = 13;

    // Create a nonsense object just to check the complementarity.
    const vector<vector<size_t>> owned_cells_per_rank {
        vector<size_t> { 0, 5, 7 },
        vector<size_t> { 2, 3 }
    };

    const auto non_owned_cells_per_rank = get_non_owned_cells_per_rank(
        owned_cells_per_rank, num_cells
    );

    for (size_t rank = 0; rank < 2; ++rank)
    {
        const auto& owned_cells = owned_cells_per_rank.at(rank);
        const auto& non_owned = non_owned_cells_per_rank.at(rank);

        ASSERT_EQ(owned_cells.size() + non_owned.size(), num_cells,
            "the owned and non-owned cells are not the same "
            "total number as the existing cells");

        for (size_t i = 0; i < num_cells; ++i)
        {
            if (in_set(i, owned_cells))
            {
                ASSERT(!in_set(i, non_owned),
                    "a cell is in both the owned and non-owned lists");
            }
            else
            {
                ASSERT(in_set(i, non_owned),
                    "a cell is not in either the owned or non-owned lists");
            }
        }
    }
)

ADD_TEST(test_get_cell_ownership_per_rank_metadata,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    const size_t num_ranks = 4;
    const size_t num_cells = 12;

    const auto expected_owned_cells_per_rank
        = get_owned_cells_per_rank(num_cells, num_ranks);
    const auto expected_cell_parent_mpi_ranks
        = get_cell_mpi_ranks(expected_owned_cells_per_rank);
    const auto expected_non_owned_cells_per_rank = get_non_owned_cells_per_rank(
        expected_owned_cells_per_rank, num_cells
    );

    fill_mpi_rank_and_cell_ownership(mpi_comm, num_cells);

    for (size_t rank = 0; rank < num_ranks; ++rank)
    {
        ASSERT_EQ_VEC(
            mpi_comm.mpi_rank_owned_cells.at(rank),
            expected_owned_cells_per_rank.at(rank),
            "the right owned cells per ranks were not filled in for all ranks"
        );

        ASSERT_EQ_VEC(
            mpi_comm.mpi_rank_non_owned_cells.at(rank),
            expected_non_owned_cells_per_rank.at(rank),
            "the right non-owned cells per ranks were not filled in "
            "for all ranks"
        );
    }

    ASSERT_EQ_VEC(
        mpi_comm.cell_parent_mpi_ranks, expected_cell_parent_mpi_ranks,
        "the right cell parents were not filled in"
    );

)

ADD_TEST(test_create_all_cells_for_all_mpi_ranks_with_correct_fields,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    System system;

    const RVec cell_size { 1.0, 2.0, 3.0 };
    const RVec origin0 { 0.1, 0.2, 0.3 };
    const RVec origin1 { 0.4, 0.5, 0.6 };
    const RVec origin2 { 0.7, 0.8, 0.9 };

    const vector<size_t> list0_neighbours {1, 2, 3};
    const vector<size_t> list1_neighbours {0, 2};
    const vector<size_t> list2_neighbours;

    if (is_master(mpi_comm))
    {
        CellList list0 { 0, origin0, cell_size };
        CellList list1 { 0, origin1, cell_size };
        CellList list2 { 0, origin2, cell_size };

        list0.to_neighbours = list0_neighbours;
        list1.to_neighbours = list1_neighbours;
        list2.to_neighbours = list2_neighbours;

        system.cell_lists.push_back(list0);
        system.cell_lists.push_back(list1);
        system.cell_lists.push_back(list2);
    }

    mpi_create_all_cell_lists(system, 3, mpi_comm);

    ASSERT_EQ(system.cell_lists.size(), 3,
        "not all ranks have the correct number of cells");

    ASSERT_EQ_VEC(system.cell_lists[0].origin, origin0,
        "not all ranks have the correct origin of cell 0");
    ASSERT_EQ_VEC(system.cell_lists[1].origin, origin1,
        "not all ranks have the correct origin of cell 1");
    ASSERT_EQ_VEC(system.cell_lists[2].origin, origin2,
        "not all ranks have the correct origin of cell 2");

    ASSERT_EQ_VEC(system.cell_lists[0].to_neighbours, list0_neighbours,
        "not all ranks have the correct neighbours list of cell 0");
    ASSERT_EQ_VEC(system.cell_lists[1].to_neighbours, list1_neighbours,
        "not all ranks have the correct neighbours list of cell 1");
    ASSERT_EQ_VEC(system.cell_lists[2].to_neighbours, list2_neighbours,
        "not all ranks have the correct neighbours list of cell 2");

    for (const auto list : system.cell_lists)
    {
        ASSERT_EQ_VEC(list.size, cell_size,
            "not all ranks and cells have the correct cell size");
    }
)

ADD_TEST(test_divide_cell_lists_onto_proper_ranks,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    System system;

    const vector<real> list1_vels { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
    const vector<real> list2_vels { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
    const vector<real> list3_vels { 6.0, 7.0, 8.0 };
    const vector<real> list4_vels { 6.1, 6.2, 6.3 };

    if (is_master(mpi_comm))
    {
        CellList list1 {2, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};
        CellList list2 {3, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};
        CellList list3 {1, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};
        CellList list4 {1, RVec {0.0, 0.0, 0.0}, RVec {0.0, 0.0, 0.0}};

        list1.add_atom(0.0, 0.0, 0.0);
        list1.add_atom(1.0, 0.0, 0.0);
        ASSERT_EQ(list1.xs.size(), list1_vels.size(),
            "bad setup: the position and velocity vectors are not of "
            "same size for list 1");
        list1.vs = list1_vels;

        list2.add_atom(0.0, 0.0, 0.0);
        list2.add_atom(1.0, 0.0, 0.0);
        list2.add_atom(2.0, 0.0, 0.0);
        ASSERT_EQ(list2.xs.size(), list2_vels.size(),
            "bad setup: the position and velocity vectors are not of "
            "same size for list 2");
        list2.vs = list2_vels;

        list3.add_atom(0.0, 1.0, 0.0);
        ASSERT_EQ(list3.xs.size(), list3_vels.size(),
            "bad setup: the position and velocity vectors are not of "
            "same size for list 3");
        list3.vs = list3_vels;

        list4.add_atom(0.0, 0.0, 1.0);
        ASSERT_EQ(list4.xs.size(), list4_vels.size(),
            "bad setup: the position and velocity vectors are not of "
            "same size for list 4");
        list4.vs = list4_vels;

        system.cell_lists.push_back(list1);
        system.cell_lists.push_back(list2);
        system.cell_lists.push_back(list3);
        system.cell_lists.push_back(list4);
    }

    mpi_init_cell_lists_and_transfer(system, mpi_comm);

    // Check that each rank has the number of atoms of its cell
    switch (mpi_comm.rank)
    {
        case 0:
            ASSERT_EQ(system.num_atoms(), 2,
                "rank 0 does not have the correct number of atoms");
            break;

        case 1:
            ASSERT_EQ(system.num_atoms(), 3,
                "rank 1 does not have the correct number of atoms");
            break;

        case 2:
            ASSERT_EQ(system.num_atoms(), 1,
                "rank 2 does not have the correct number of atoms");
            break;

        case 3:
            ASSERT_EQ(system.num_atoms(), 1,
                "rank 3 does not have the correct number of atoms");
            break;
        default:
            break;
    }

    // Now check that all ranks got the correct atom positions and velocities,
    // and that the forces are initialized to 0
    vector<real> xs, fs;
    switch (mpi_comm.rank)
    {
        case 0:
            xs = vector<real> { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
            fs = vector<real> (xs.size(), 0.0);

            ASSERT_EQ(system.cell_lists[0].num_atoms(), 2,
                "rank 0 does not have the correct number of atoms");
            ASSERT_EQ_VEC(system.cell_lists[0].xs, xs,
                "rank 0 did not keep its cell positions correctly");
            ASSERT_EQ_VEC(system.cell_lists[0].vs, list1_vels,
                "rank 0 did not keep its cell velocities correctly");
            ASSERT_EQ_VEC(system.cell_lists[0].fs, fs,
                "rank 0 did not keep the correct number of zero-initialized forces");
            break;

        case 1:
            xs = vector<real> { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0 };
            fs = vector<real> (xs.size(), 0.0);

            ASSERT_EQ(system.cell_lists[1].num_atoms(), 3,
                "rank 1 does not have the correct number of atoms");
            ASSERT_EQ_VEC(system.cell_lists[1].xs, xs,
                "rank 1 did not get its positions correctly");
            ASSERT_EQ_VEC(system.cell_lists[1].vs, list2_vels,
                "rank 1 did not get its velocities correctly");
            ASSERT_EQ_VEC(system.cell_lists[1].fs, fs,
                "rank 1 did not get the correct number of zero-initialized forces");
            break;

        case 2:
            xs = vector<real> { 0.0, 1.0, 0.0 };
            fs = vector<real> (xs.size(), 0.0);

            ASSERT_EQ(system.cell_lists[2].num_atoms(), 1,
                "rank 2 does not have the correct number of atoms");
            ASSERT_EQ_VEC(system.cell_lists[2].xs, xs,
                "rank 2 did not get its positions correctly");
            ASSERT_EQ_VEC(system.cell_lists[2].vs, list3_vels,
                "rank 2 did not get its velocities correctly");
            ASSERT_EQ_VEC(system.cell_lists[2].fs, fs,
                "rank 2 did not get the correct number of zero-initialized forces");
            break;

        case 3:
            xs = vector<real> { 0.0, 0.0, 1.0 };
            fs = vector<real> (xs.size(), 0.0);

            ASSERT_EQ(system.cell_lists[3].num_atoms(), 1,
                "rank 3 does not have the correct number of atoms");
            ASSERT_EQ_VEC(system.cell_lists[3].xs, xs,
                "rank 3 did not get its positions correctly");
            ASSERT_EQ_VEC(system.cell_lists[3].vs, list4_vels,
                "rank 3 did not get its velocities correctly");
            ASSERT_EQ_VEC(system.cell_lists[3].fs, fs,
                "rank 3 did not get the correct number of zero-initialized forces");
            break;
    }
)

ADD_TEST(test_init_mpi_for_test_suite,
    int num_ranks;

    ASSERT_EQ(MPI_Comm_size(MPI_COMM_WORLD, &num_ranks), MPI_SUCCESS,
        "MPI could not be initialized for the testing suite");

    ASSERT_EQ(num_ranks, 4, "the number of MPI ranks for testing has to be 4");
)

ADD_TEST(test_create_cell_mpi_communication_groups_for_sending,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // 4 ranks with a few cells each
    mpi_comm.cell_parent_mpi_ranks = vector<size_t> { 3, 2, 1, 0, 3, 2 };

    System system;

    CellList list1 {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list2 {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list3 {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list4 {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list5 {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list6 {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};

    // Neighbouring cell indices
    list1.to_neighbours = vector<size_t> { 1, 2 }; // rank 3 -> 2, 1
    list2.to_neighbours = vector<size_t> { 2, 3 }; // rank 2 -> 1, 0
    list3.to_neighbours = vector<size_t> { 3, 4 }; // rank 1 -> 0, 3
    list4.to_neighbours = vector<size_t> { 4, 5 }; // rank 0 -> 3, 2
    list5.to_neighbours = vector<size_t> { 5, 0 }; // rank 3 -> 2, 3
    list6.to_neighbours = vector<size_t> { 0, 1 }; // rank 2 -> 3, 2

    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);
    system.cell_lists.push_back(list3);
    system.cell_lists.push_back(list4);
    system.cell_lists.push_back(list5);
    system.cell_lists.push_back(list6);

    // The comm group will contain the ranks of the neighbouring cells
    // and the parent cell
    const vector<int> ranks_cell1 { 3, 2, 1 };
    const vector<int> ranks_cell2 { 2, 1, 0 };
    const vector<int> ranks_cell3 { 1, 0, 3 };
    const vector<int> ranks_cell4 { 0, 3, 2 };

    // These two cells have double members (their neighbours include their
    // parent rank), which are removed
    const vector<int> ranks_cell5 { 3, 2 };
    const vector<int> ranks_cell6 { 2, 3 };

    MPI_Group expected_group1, expected_group2, expected_group3,
              expected_group4, expected_group5, expected_group6,
              group1, group2, group3, group4, group5, group6,
              world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    // To test the results, we create the expected groups manually
    MPI_Group_incl(world_group, ranks_cell1.size(), ranks_cell1.data(),
        &expected_group1);
    MPI_Group_incl(world_group, ranks_cell2.size(), ranks_cell2.data(),
        &expected_group2);
    MPI_Group_incl(world_group, ranks_cell3.size(), ranks_cell3.data(),
        &expected_group3);
    MPI_Group_incl(world_group, ranks_cell4.size(), ranks_cell4.data(),
        &expected_group4);
    MPI_Group_incl(world_group, ranks_cell5.size(), ranks_cell5.data(),
        &expected_group5);
    MPI_Group_incl(world_group, ranks_cell6.size(), ranks_cell6.data(),
        &expected_group6);

    // Now extract them from the MPI_Comm records in the result and compare
    int result = MPI_UNEQUAL, cell_comm_root = -1;
    mpi_create_cell_comm_groups(mpi_comm, system.cell_lists);
    const auto& comm_groups = mpi_comm.cell_comm_groups;

    // Use the (known) ranks for each cell to check the groups
    // Owner: 3
    if (is_rank(3, mpi_comm) || is_rank(2, mpi_comm) || is_rank(1, mpi_comm))
    {
        // Get the group from the created MPI_Comm record
        MPI_Comm_group(comm_groups[0].comm, &group1);

        // Compare the group to the expected result.
        MPI_Group_compare(group1, expected_group1, &result);

        // Ensure that they are not unequal.
        ASSERT(static_cast<bool>(result != MPI_UNEQUAL),
               "cell 0 does not have the correct ranks in the comm group");

        // Ensure that the sending ranks match those recorded.
        const set<size_t> to_ranks_expected { 2, 1 };
        const set<size_t> to_ranks(
            comm_groups[0].to_ranks.cbegin(), comm_groups[0].to_ranks.cend()
        );

        ASSERT_EQ_VEC(to_ranks, to_ranks_expected,
                      "cell 0 does not have the correct to_ranks list");

        // Check that the root corresponds to the owning rank in MPI_COMM_WORLD:
        // get the number by translating from MPI_COMM_WORLD to the group comm.
        const int owner = 3;
        MPI_Group_translate_ranks(world_group, 1, &owner,
            group1, &cell_comm_root);

        ASSERT_EQ(cell_comm_root, comm_groups[0].root,
            "cell 0 does not have the correct root handle");
    }

    // Owner: 2
    if (is_rank(2, mpi_comm) || is_rank(1, mpi_comm) || is_rank(0, mpi_comm))
    {
        MPI_Comm_group(comm_groups[1].comm, &group2);
        MPI_Group_compare(group2, expected_group2, &result);
        ASSERT(static_cast<bool>(result != MPI_UNEQUAL),
               "cell 1 does not have the correct ranks in the comm group");

        const set<size_t> to_ranks_expected { 1, 0 };
        const set<size_t> to_ranks(
            comm_groups[1].to_ranks.cbegin(), comm_groups[1].to_ranks.cend()
        );
        ASSERT_EQ_VEC(to_ranks, to_ranks_expected,
                      "cell 1 does not have the correct to_ranks list");

        const int owner = 2;
        MPI_Group_translate_ranks(world_group, 1, &owner,
            group2, &cell_comm_root);

        ASSERT_EQ(cell_comm_root, comm_groups[1].root,
            "cell 1 does not have the correct root handle");
    }

    // Owner: 1
    if (is_rank(1, mpi_comm) || is_rank(0, mpi_comm) || is_rank(3, mpi_comm))
    {
        MPI_Comm_group(comm_groups[2].comm, &group3);
        MPI_Group_compare(group3, expected_group3, &result);
        ASSERT(static_cast<bool>(result != MPI_UNEQUAL),
               "cell 2 does not have the correct ranks in the comm group");

        // Ensure that the sending ranks match those recorded.
        const set<size_t> to_ranks_expected { 0, 3 };
        const set<size_t> to_ranks(
            comm_groups[2].to_ranks.cbegin(), comm_groups[2].to_ranks.cend()
        );
        ASSERT_EQ_VEC(to_ranks, to_ranks_expected,
                      "cell 2 does not have the correct to_ranks list");

        const int owner = 1;
        MPI_Group_translate_ranks(world_group, 1, &owner,
            group3, &cell_comm_root);

        ASSERT_EQ(cell_comm_root, comm_groups[2].root,
            "cell 2 does not have the correct root handle");
    }

    // Owner: 0
    if (is_rank(0, mpi_comm) || is_rank(3, mpi_comm) || is_rank(2, mpi_comm))
    {
        MPI_Comm_group(comm_groups[3].comm, &group4);
        MPI_Group_compare(group4, expected_group4, &result);
        ASSERT(static_cast<bool>(result != MPI_UNEQUAL),
               "cell 3 does not have the correct ranks in the comm group");

        const set<size_t> to_ranks_expected { 3, 2 };
        const set<size_t> to_ranks(
            comm_groups[3].to_ranks.cbegin(), comm_groups[3].to_ranks.cend()
        );
        ASSERT_EQ_VEC(to_ranks, to_ranks_expected,
                      "cell 3 does not have the correct to_ranks list");

        const int owner = 0;
        MPI_Group_translate_ranks(world_group, 1, &owner,
            group4, &cell_comm_root);

        ASSERT_EQ(cell_comm_root, comm_groups[3].root,
            "cell 3 does not have the correct root handle");
    }

    // Owner: 3
    if (is_rank(3, mpi_comm) || is_rank(2, mpi_comm))
    {
        MPI_Comm_group(comm_groups[4].comm, &group5);
        MPI_Group_compare(group5, expected_group5, &result);
        ASSERT(static_cast<bool>(result != MPI_UNEQUAL),
               "cell 4 does not have the correct ranks in the comm group");

        const set<size_t> to_ranks_expected { 2 };
        const set<size_t> to_ranks(
            comm_groups[4].to_ranks.cbegin(), comm_groups[4].to_ranks.cend()
        );
        ASSERT_EQ_VEC(to_ranks, to_ranks_expected,
                      "cell 4 does not have the correct to_ranks list");

        const int owner = 3;
        MPI_Group_translate_ranks(world_group, 1, &owner,
            group5, &cell_comm_root);

        ASSERT_EQ(cell_comm_root, comm_groups[4].root,
            "cell 4 does not have the correct root handle");
    }

    // Owner: 2
    if (is_rank(2, mpi_comm) || is_rank(3, mpi_comm))
    {
        MPI_Comm_group(comm_groups[5].comm, &group6);
        MPI_Group_compare(group6, expected_group6, &result);
        ASSERT(static_cast<bool>(result != MPI_UNEQUAL),
               "cell 5 does not have the correct ranks in the comm group");

        const set<size_t> to_ranks_expected { 3 };
        const set<size_t> to_ranks(
            comm_groups[5].to_ranks.cbegin(), comm_groups[5].to_ranks.cend()
        );
        ASSERT_EQ_VEC(to_ranks, to_ranks_expected,
                      "cell 5 does not have the correct to_ranks list");

        const int owner = 2;
        MPI_Group_translate_ranks(world_group, 1, &owner,
            group6, &cell_comm_root);

        ASSERT_EQ(cell_comm_root, comm_groups[5].root,
            "cell 5 does not have the correct root handle");
    }

    MPI_Group_free(&world_group);
)

ADD_TEST(test_create_record_of_which_cells_every_rank_will_receive,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // 4 ranks with a few cells each
    mpi_comm.cell_parent_mpi_ranks = vector<size_t> { 3, 2, 1, 0, 3, 2 };

    System system;

    CellList list0 {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list1 {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list2 {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list3 {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list4 {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list5 {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};

    // Neighbouring cell indices
    list0.to_neighbours = vector<size_t> { 1, 2 }; // rank 3 -> 2, 1
    list1.to_neighbours = vector<size_t> { 2, 3 }; // rank 2 -> 1, 0
    list2.to_neighbours = vector<size_t> { 3, 4 }; // rank 1 -> 0, 3
    list3.to_neighbours = vector<size_t> { 4, 5 }; // rank 0 -> 3, 2
    list4.to_neighbours = vector<size_t> { 5, 0 }; // rank 3 -> 2, 3
    list5.to_neighbours = vector<size_t> { 0, 1 }; // rank 2 -> 3, 2

    system.cell_lists.push_back(list0);
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);
    system.cell_lists.push_back(list3);
    system.cell_lists.push_back(list4);
    system.cell_lists.push_back(list5);

    // The ranks will be receiving these cells, based on the to_neighbours list
    const set<size_t> rank0_expects { 1, 2 };
    const set<size_t> rank1_expects { 0, 1 };

    // cell no. 5 removed, since it is owned by the rank
    const set<size_t> rank2_expects { 0, 3, 4 };

    // cell no. 4 removed, since it is owned by the rank
    const set<size_t> rank3_expects { 2, 3, 5 };

    mpi_create_sending_and_receiving_cell_lists_for_ranks(
        mpi_comm, system.cell_lists
    );

    const auto& received = mpi_comm.mpi_rank_received_cells;
    ASSERT_EQ(received.size(), mpi_comm.num_ranks,
        "mismatch between number of ranks and number of received cells per rank");

    // Use sets since we don't care about the order
    const set<size_t> rank0_results(received[0].cbegin(), received[0].cend());
    const set<size_t> rank1_results(received[1].cbegin(), received[1].cend());
    const set<size_t> rank2_results(received[2].cbegin(), received[2].cend());
    const set<size_t> rank3_results(received[3].cbegin(), received[3].cend());

    ASSERT_EQ_VEC(rank0_results, rank0_expects,
        "rank 0 does not have the correct received cells");
    ASSERT_EQ_VEC(rank1_results, rank1_expects,
        "rank 1 does not have the correct received cells");
    ASSERT_EQ_VEC(rank2_results, rank2_expects,
        "rank 2 does not have the correct received cells");
    ASSERT_EQ_VEC(rank3_results, rank3_expects,
        "rank 3 does not have the correct received cells");
)

ADD_TEST(test_create_record_of_which_cells_every_rank_will_send,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // 4 ranks with a few cells each
    mpi_comm.cell_parent_mpi_ranks = vector<size_t> { 0, 0, 0, 1, 1, 1 };

    System system;

    CellList list0 {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list1 {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list2 {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list3 {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list4 {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list5 {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};

    // Neighbouring cell indices
    list0.to_neighbours = vector<size_t> { 1, 2 }; // rank 0 -> 0, 0 (no send)
    list1.to_neighbours = vector<size_t> { 2, 3 }; // rank 0 -> 0, 1 (send)
    list2.to_neighbours = vector<size_t> { 3, 4 }; // rank 0 -> 1, 1 (send)
    list3.to_neighbours = vector<size_t> { 4, 5 }; // rank 1 -> 1, 1 (no send)
    list4.to_neighbours = vector<size_t> { 5, 0 }; // rank 1 -> 1, 0 (send)
    list5.to_neighbours = vector<size_t> { 0, 1 }; // rank 1 -> 0, 0 (send)

    system.cell_lists.push_back(list0);
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);
    system.cell_lists.push_back(list3);
    system.cell_lists.push_back(list4);
    system.cell_lists.push_back(list5);

    // The ranks will be sending these cells, based on the to_neighbours list
    const set<size_t> rank0_expects { 1, 2 };
    const set<size_t> rank1_expects { 4, 5 };

    // These ranks will not be sending any cells
    const set<size_t> rank2_expects { };
    const set<size_t> rank3_expects { };

    // mpi_create_sending_cell_lists_for_ranks(mpi_comm, system.cell_lists);
    mpi_create_sending_and_receiving_cell_lists_for_ranks(
        mpi_comm, system.cell_lists
    );

    const auto& sending = mpi_comm.mpi_rank_sending_cells;
    ASSERT_EQ(sending.size(), mpi_comm.num_ranks,
        "the sending list does not contain one list per MPI rank");

    // Use sets since we don't care about the order
    const set<size_t> rank0_results(sending[0].cbegin(), sending[0].cend());
    const set<size_t> rank1_results(sending[1].cbegin(), sending[1].cend());
    const set<size_t> rank2_results(sending[2].cbegin(), sending[2].cend());
    const set<size_t> rank3_results(sending[3].cbegin(), sending[3].cend());

    ASSERT_EQ_VEC(rank0_results, rank0_expects,
        "rank 0 does not have the correct sending cells");
    ASSERT_EQ_VEC(rank1_results, rank1_expects,
        "rank 1 does not have the correct sending cells");
    ASSERT_EQ_VEC(rank2_results, rank2_expects,
        "rank 2 does not have the correct sending cells");
    ASSERT_EQ_VEC(rank3_results, rank3_expects,
        "rank 3 does not have the correct sending cells");
)

ADD_TEST(test_sending_and_receiving_cell_lists_are_consistent,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // 4 ranks with a few cells each
    mpi_comm.cell_parent_mpi_ranks = vector<size_t> { 3, 2, 1, 0, 3, 2 };

    System system;

    CellList list0 {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list1 {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list2 {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list3 {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list4 {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list5 {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};

    // Neighbouring cell indices
    list0.to_neighbours = vector<size_t> { 1, 2 }; // rank 3 -> 2, 1
    list1.to_neighbours = vector<size_t> { 2, 3 }; // rank 2 -> 1, 0
    list2.to_neighbours = vector<size_t> { 3, 4 }; // rank 1 -> 0, 3
    list3.to_neighbours = vector<size_t> { 4, 5 }; // rank 0 -> 3, 2
    list4.to_neighbours = vector<size_t> { 5, 0 }; // rank 3 -> 2, 3
    list5.to_neighbours = vector<size_t> { 0, 1 }; // rank 2 -> 3, 2

    // And their corresponding parent ranks
    vector<vector<size_t>> to_neighbour_ranks {
        vector<size_t> { 2, 1 },
        vector<size_t> { 1, 0 },
        vector<size_t> { 0, 3 },
        vector<size_t> { 3, 2 },
        vector<size_t> { 2 },
        vector<size_t> { 3 }
    };

    system.cell_lists.push_back(list0);
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);
    system.cell_lists.push_back(list3);
    system.cell_lists.push_back(list4);
    system.cell_lists.push_back(list5);

    mpi_create_sending_and_receiving_cell_lists_for_ranks(
        mpi_comm, system.cell_lists
    );

    const auto& received = mpi_comm.mpi_rank_received_cells;
    const auto& sending = mpi_comm.mpi_rank_sending_cells;

    const vector<set<size_t>> recv_sets {
        set<size_t> (received[0].cbegin(), received[0].cend()),
        set<size_t> (received[1].cbegin(), received[1].cend()),
        set<size_t> (received[2].cbegin(), received[2].cend()),
        set<size_t> (received[3].cbegin(), received[3].cend())
    };

    const vector<set<size_t>> send_sets {
        set<size_t> (sending[0].cbegin(), sending[0].cend()),
        set<size_t> (sending[1].cbegin(), sending[1].cend()),
        set<size_t> (sending[2].cbegin(), sending[2].cend()),
        set<size_t> (sending[3].cbegin(), sending[3].cend())
    };

    // For every cell that is send from every rank, ensure that it is received
    // by a corresponding one.
    for (const auto sending_cells : send_sets)
    {
        for (const auto i : sending_cells)
        {
            const auto from_rank = mpi_comm.cell_parent_mpi_ranks.at(i);

            for (const auto to_rank : to_neighbour_ranks.at(i))
            {
                ASSERT(static_cast<bool>(from_rank != to_rank),
                       "a cell will be sent to itself");

                ASSERT(in_set(i, recv_sets.at(to_rank)),
                       "a sent cell was not found in the receiver");
            }
        }
    }

    // And vice versa: for every cell that is received, from every rank,
    // ensure that it is sent by the corresponding one.
    for (auto to_rank = 0; to_rank < recv_sets.size(); ++to_rank)
    {
        const auto& receiving_cells = recv_sets.at(to_rank);

        for (const auto i : receiving_cells)
        {
            const auto from_rank = mpi_comm.cell_parent_mpi_ranks.at(i);

            ASSERT(static_cast<bool>(from_rank != to_rank),
                   "a cell will be received by itself");

            ASSERT(in_set(i, send_sets.at(from_rank)),
                   "a received cell was not found in the sender");
        }
    }
)

ADD_TEST(test_fill_in_communicators_and_sendrecv_cell_lists,
    MPIRank mpi_comm, mpi_comm_expected;
    init_MPI(mpi_comm);
    init_MPI(mpi_comm_expected);

    mpi_comm.cell_parent_mpi_ranks = vector<size_t> { 3, 2, 1, 0, 3, 2 };
    mpi_comm_expected.cell_parent_mpi_ranks = vector<size_t> (
        mpi_comm.cell_parent_mpi_ranks.cbegin(),
        mpi_comm.cell_parent_mpi_ranks.cend()
    );

    System system;

    system.cell_lists = {
        CellList {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}}
    };

    // Neighbouring cell indices
    system.cell_lists[0].to_neighbours = vector<size_t> { 1, 2 };
    system.cell_lists[1].to_neighbours = vector<size_t> { 2, 3 };
    system.cell_lists[2].to_neighbours = vector<size_t> { 3, 4 };
    system.cell_lists[3].to_neighbours = vector<size_t> { 4, 5 };
    system.cell_lists[4].to_neighbours = vector<size_t> { 5, 0 };
    system.cell_lists[5].to_neighbours = vector<size_t> { 0, 1 };

    // Constuct comparison object by running the expected functions
    // one-by-one
    mpi_create_cell_comm_groups(mpi_comm_expected, system.cell_lists);
    mpi_create_sending_and_receiving_cell_lists_for_ranks(
        mpi_comm_expected, system.cell_lists
    );

    mpi_fill_communication_data(mpi_comm, system);

    ASSERT_EQ(mpi_comm.cell_comm_groups.size(),
              mpi_comm_expected.cell_comm_groups.size(),
              "the comm groups were not filled in");

    for (size_t i = 0; i < system.cell_lists.size(); ++i)
    {
        const auto& comm_group = mpi_comm.cell_comm_groups.at(i);
        const auto& comm_group_expected
            = mpi_comm_expected.cell_comm_groups.at(i);

        ASSERT_EQ(comm_group.root, comm_group_expected.root,
            "root in a cell MPICellComm object was not matching");
        ASSERT_EQ_VEC(comm_group.to_ranks, comm_group_expected.to_ranks,
            "to_ranks in a cell MPICellComm object was not matching");

        if (comm_group.comm != MPI_COMM_NULL && comm_group_expected.comm != MPI_COMM_NULL)
        {
            int result;
            const auto comm_groups_equal = static_cast<bool>(
                MPI_Comm_compare(comm_group.comm, comm_group_expected.comm, &result)
            );
            ASSERT(static_cast<bool>(result != MPI_UNEQUAL),
                "comm in a cell MPICellComm object was not matching");
        }
    }

    for (size_t rank = 0; rank < 4; ++rank)
    {
        const auto& recv_cells = mpi_comm.mpi_rank_received_cells.at(rank);
        const auto& send_cells = mpi_comm.mpi_rank_sending_cells.at(rank);

        const auto& recv_cells_expected
            = mpi_comm_expected.mpi_rank_received_cells.at(rank);
        const auto& send_cells_expected
            = mpi_comm_expected.mpi_rank_sending_cells.at(rank);

        ASSERT_EQ_VEC(recv_cells, recv_cells_expected,
            "receiving cells were not matching for a rank");
        ASSERT_EQ_VEC(send_cells, send_cells_expected,
            "receiving cells were not matching for a rank");
    }
)

template <typename T>
static size_t find_value_index(const T needle, const vector<T>& haystack)
{
    const auto it = find(haystack.cbegin(), haystack.cend(), needle);
    return distance(haystack.cbegin(), it);
}

ADD_TEST(test_sync_num_of_transmitted_positions,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // 4 ranks with a few cells each
    mpi_comm.cell_parent_mpi_ranks = vector<size_t> { 3, 2, 1, 0, 3, 2 };

    System system;

    // Create a couple of atoms in the first list, owned by rank 3
    CellList list0 {2, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list1 {1, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list2 {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list3 {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list4 {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list5 {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};

    // Neighbouring cell indices:
    // Two atoms in list0 will be sent from rank 3 to ranks 2 and 1
    list0.to_neighbours = vector<size_t> { 1, 2 }; // rank 3 -> 2, 1

    // One atom in list1 will be sent from rank 2 to ranks 1 and 0
    list1.to_neighbours = vector<size_t> { 2, 3 };

    // Empty
    list2.to_neighbours = vector<size_t> { 3, 4 };
    list3.to_neighbours = vector<size_t> { 4, 5 };
    list4.to_neighbours = vector<size_t> { 5, 0 };
    list5.to_neighbours = vector<size_t> { 0, 1 };

    if (is_rank(3, mpi_comm))
    {
        list0.add_atom(0.0, 1.0, 2.0);
        list0.add_atom(3.0, 4.0, 5.0);

        ASSERT_EQ(list0.num_atoms(), 2,
            "rank 3 cell 0 does not start with the correct number of atoms");
    }
    else
    {
        ASSERT_EQ(list0.xs.size(), 0,
            "other ranks do not start with an empty cell0");
    }

    if (is_rank(2, mpi_comm))
    {
        list1.add_atom(0.0, 1.0, 2.0);

        ASSERT_EQ(list1.num_atoms(), 1,
            "rank 2 cell 1 does not start with the correct number of atoms");
    }
    else
    {
        ASSERT_EQ(list1.xs.size(), 0,
            "other ranks do not start with an empty cell1");
    }

    system.cell_lists.push_back(list0);
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);
    system.cell_lists.push_back(list3);
    system.cell_lists.push_back(list4);
    system.cell_lists.push_back(list5);


    // Setup the communication records
    mpi_create_cell_comm_groups(mpi_comm, system.cell_lists);
    mpi_create_sending_and_receiving_cell_lists_for_ranks(
        mpi_comm, system.cell_lists
    );

    // Get the number of positions to receive per cell
    const auto num_recv_per_cell = mpi_sync_number_of_transmitted_atoms(
        system.cell_lists, mpi_comm
    );

    // Assert that we have the correct number for every rank
    // Rank 0 receives cell1 from rank 2
    if (is_rank(0, mpi_comm))
    {
        // In the mpi_rank_received_cells for this rank, find the index
        // corresponding to cell1. Then check that position in the returned
        // num_recv_per_cell vector for correctness, since they should be
        // identically indexed.

        // First, assert that cell1 is in the received list for this rank.
        const auto& recv_cells = mpi_comm.mpi_rank_received_cells.at(0);
        const auto result = in_set(static_cast<size_t>(1), recv_cells);
        ASSERT(result, "rank 0: expected to receive cell1, but did not");

        if (result)
        {
            // Then, get its index.
            const auto icell1 = find_value_index(
                static_cast<size_t>(1), recv_cells
            );

            // Loop through the indices in the vector and for the matching
            // index, check the number. Assert that the others are empty.
            for (int j = 0; j < recv_cells.size(); ++j)
            {
                const auto& num = num_recv_per_cell.at(j);

                if (j == icell1)
                {
                    ASSERT_EQ(
                        num, 1,
                        "rank 0 receives the wrong number of cell1 from "
                        "rank 2"
                    );
                }
                else
                {
                    ASSERT_EQ(num, 0, "rank 0 will receive unexpected cells");
                }
            }
        }
    }

    // Rank 1 receives cell0 from rank 3 and cell1 from rank 2
    if (is_rank(1, mpi_comm))
    {
        // Two received cells, ensure that both are present.
        const auto& recv_cells = mpi_comm.mpi_rank_received_cells.at(1);
        const auto result0 = in_set(static_cast<size_t>(0), recv_cells);
        const auto result1 = in_set(static_cast<size_t>(1), recv_cells);

        ASSERT(result0 && result1,
            "rank 1: expected to receive cell0 and cell1, but did not");

        if (result0 && result1)
        {
            // Then, get both indices.
            const auto icell0 = find_value_index(
                static_cast<size_t>(0), recv_cells
            );
            const auto icell1 = find_value_index(
                static_cast<size_t>(1), recv_cells
            );

            // Loop through the indices in the vector and for the matching
            // index, check the number. Assert that the others are empty.
            for (int j = 0; j < recv_cells.size(); ++j)
            {
                const auto& num = num_recv_per_cell.at(j);

                if (j == icell0)
                {
                    ASSERT_EQ(
                        num, 2,
                        "rank 1 receives the wrong number of cell0 from "
                        "rank 3"
                    );
                }
                else if (j == icell1)
                {
                    ASSERT_EQ(
                        num, 1,
                        "rank 1 receives the wrong number of cell1 from "
                        "rank 2"
                    );
                }
                else
                {
                    char buf[256];
                    snprintf(buf, 256,
                        "rank 1 will receive %d atoms in cell %d, "
                        "but expected none", num, recv_cells.at(j));
                    ASSERT_EQ(num, 0, buf);
                }
            }
        }
    }

    // Rank 2 receives cell0 from rank 3
    if (is_rank(2, mpi_comm))
    {
        const auto& recv_cells = mpi_comm.mpi_rank_received_cells.at(2);
        const auto result = in_set(static_cast<size_t>(0), recv_cells);
        ASSERT(result, "rank 2: expected to receive cell0, but did not");

        if (result)
        {
            // Then, get its index.
            const auto icell0 = find_value_index(
                static_cast<size_t>(0), recv_cells
            );

            // Loop through the indices in the vector and for the matching
            // index, check the number. Assert that the others are empty.
            for (int j = 0; j < recv_cells.size(); ++j)
            {
                const auto& num = num_recv_per_cell.at(j);

                if (j == icell0)
                {
                    ASSERT_EQ(
                        num, 2,
                        "rank 2 receives the wrong number of cell0 from "
                        "rank 3"
                    );
                }
                else
                {
                    ASSERT_EQ(num, 0, "rank 0 will receive unexpected cells");
                }
            }
        }
    }

    // Rank 3 does not receive any cells
    if (is_rank(3, mpi_comm))
    {
        const auto& recv_cells = mpi_comm.mpi_rank_received_cells.at(3);
        ASSERT_EQ(num_recv_per_cell.size(), recv_cells.size(),
            "rank 3: incorrect size of returned array");

        for (const auto num : num_recv_per_cell)
        {
            ASSERT_EQ(num, 0, "rank 3 will receive unexpected cells");
        }
    }
)

ADD_TEST(test_reserve_memory_for_recv_cell_lists_positions_and_forces,
    // Corresponds to the field in MPIRank
    const std::vector<size_t> mpi_rank_received_cells { 3, 2, 5 };

    // Cell 3: 5 atoms
    // Cell 2: 8 atoms
    // Cell 5: 1 atom
    // Others are empty
    const std::vector<size_t> num_recv_per_cell { 5, 8, 1 };

    vector<CellList> cell_lists {
        CellList {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}}
    };

    reserve_memory_for_received_lists(
        cell_lists, num_recv_per_cell, mpi_rank_received_cells
    );

    // Cells 0 and 1 receives 0 atoms
    size_t num = 0;
    ASSERT_EQ(cell_lists[0].num_atoms(), num,
              "cell0 has unexpected atoms");
    ASSERT_EQ(cell_lists[0].xs.size(), num * NDIM,
              "cell0 has unexpected positions");
    ASSERT_EQ(cell_lists[0].vs.size(), 0,
              "cell0 has unexpected velocities");
    ASSERT_EQ(cell_lists[0].fs.size(), num * NDIM,
              "cell0 has unexpected forces");

    ASSERT_EQ(cell_lists[1].num_atoms(), num,
              "cell1 has unexpected atoms");
    ASSERT_EQ(cell_lists[1].xs.size(), num * NDIM,
              "cell1 has unexpected positions");
    ASSERT_EQ(cell_lists[1].vs.size(), 0,
              "cell1 has unexpected velocities");
    ASSERT_EQ(cell_lists[1].fs.size(), num * NDIM,
              "cell1 has unexpected forces");

    // Cell 2 receives 8 atoms
    num = 8;
    ASSERT_EQ(cell_lists[2].num_atoms(), num,
              "cell2 has unexpected atoms");
    ASSERT_EQ(cell_lists[2].xs.size(), num * NDIM,
              "cell2 has unexpected positions");
    ASSERT_EQ(cell_lists[2].vs.size(), 0,
              "cell2 has unexpected velocities");
    ASSERT_EQ(cell_lists[2].fs.size(), num * NDIM,
              "cell2 has unexpected forces");

    // Cell 3 receives 5 atoms
    num = 5;
    ASSERT_EQ(cell_lists[3].num_atoms(), num,
              "cell3 has unexpected atoms");
    ASSERT_EQ(cell_lists[3].xs.size(), num * NDIM,
              "cell3 has unexpected positions");
    ASSERT_EQ(cell_lists[3].vs.size(), 0,
              "cell3 has unexpected velocities");
    ASSERT_EQ(cell_lists[3].fs.size(), num * NDIM,
              "cell3 has unexpected forces");

    // Cell 4 receives 0 atoms
    num = 0;
    ASSERT_EQ(cell_lists[4].num_atoms(), num,
              "cell4 has unexpected atoms");
    ASSERT_EQ(cell_lists[4].xs.size(), num * NDIM,
              "cell4 has unexpected positions");
    ASSERT_EQ(cell_lists[4].vs.size(), 0,
              "cell4 has unexpected velocities");
    ASSERT_EQ(cell_lists[4].fs.size(), num * NDIM,
              "cell4 has unexpected forces");

    // Cell 5 receives 1 atom
    num = 1;
    ASSERT_EQ(cell_lists[5].num_atoms(), num,
              "cell5 has unexpected atoms");
    ASSERT_EQ(cell_lists[5].xs.size(), num * NDIM,
              "cell5 has unexpected positions");
    ASSERT_EQ(cell_lists[5].vs.size(), 0,
              "cell5 has unexpected velocities");
    ASSERT_EQ(cell_lists[5].fs.size(), num * NDIM,
              "cell5 has unexpected forces");
)

ADD_TEST(test_transmitting_a_few_atom_positions_works,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // 4 ranks with a few cells each
    mpi_comm.cell_parent_mpi_ranks = vector<size_t> { 3, 2, 1, 0, 3, 2 };

    System system;

    // Create a couple of atoms in the first list, owned by rank 3
    CellList list0 {2, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list1 {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list2 {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list3 {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list4 {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list5 {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};

    // Neighbouring cell indices:
    // The atoms in list0 will be sent from rank 3 to ranks 2 and 1
    list0.to_neighbours = vector<size_t> { 1, 2 }; // rank 3 -> 2, 1
    list1.to_neighbours = vector<size_t> { 2, 3 };
    list2.to_neighbours = vector<size_t> { 3, 4 };
    list3.to_neighbours = vector<size_t> { 4, 5 };
    list4.to_neighbours = vector<size_t> { 5, 0 };
    list5.to_neighbours = vector<size_t> { 0, 1 };

    // This vector should match the positions
    const vector<real> xs_expected {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    const vector<real> fs_expected (xs_expected.size(), 0.0);

    // Add the atoms only to the owning rank
    if (is_rank(3, mpi_comm))
    {
        list0.add_atom(0.0, 1.0, 2.0);
        list0.add_atom(3.0, 4.0, 5.0);

        ASSERT_EQ_VEC(list0.xs, xs_expected,
            "the initial positions in the parent rank are incorrect");
    }
    else
    {
        ASSERT_EQ(list0.xs.size(), 0,
            "other ranks do not start with zeroed positions");
    }

    system.cell_lists.push_back(list0);
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);
    system.cell_lists.push_back(list3);
    system.cell_lists.push_back(list4);
    system.cell_lists.push_back(list5);

    // Setup the communication records
    mpi_create_cell_comm_groups(mpi_comm, system.cell_lists);
    mpi_create_sending_and_receiving_cell_lists_for_ranks(
        mpi_comm, system.cell_lists
    );

    // Transfer the atoms and verify that they (the positions) arrived
    mpi_synchronize_interaction_cell_lists(system, mpi_comm);

    // The sender should still have them
    if (is_rank(3, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ_VEC(list.xs, xs_expected,
            "the final positions in the parent rank are incorrect");
    }

    // As should rank 2
    if (is_rank(2, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ(list.num_atoms(), 2,
            "rank 2 cell does not have the correct atom count");
        ASSERT_EQ_VEC(list.xs, xs_expected,
            "rank 2 did not receive the expected positions");
        ASSERT_EQ_VEC(list.fs, fs_expected,
            "rank 2 did not created nulled-memory for the forces");
    }

    // ... and rank 1
    if (is_rank(1, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ(list.num_atoms(), 2,
            "rank 1 cell does not have the correct atom count");
        ASSERT_EQ_VEC(list.xs, xs_expected,
            "rank 1 did not receive the expected positions");
        ASSERT_EQ_VEC(list.fs, fs_expected,
            "rank 1 did not created nulled-memory for the forces");
    }

    // The remainining rank 0 should still not have any
    if (is_rank(0, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ(list0.xs.size(), 0,
            "rank 0 received positions, but it shouldn't had");
        ASSERT_EQ(list0.fs.size(), 0,
            "rank 0 received forces, but it shouldn't had");
    }

    // All the remaining cell lists should be empty, for all ranks
    for (auto i = 1; i < system.cell_lists.size(); ++i)
    {
        const auto& list = system.cell_lists.at(i);
        ASSERT_EQ(list.xs.size(), 0,
            "a rank has unexpected positions in a cell that should be empty");
        ASSERT_EQ(list.vs.size(), 0,
            "a rank has unexpected velocities in a cell that should be empty");
        ASSERT_EQ(list.fs.size(), 0,
            "a rank has unexpected forces in a cell that should be empty");
    }
)

ADD_TEST(test_transmitting_positions_erases_old_data_on_receivers,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // Same setup as the above test
    // 4 ranks with a few cells each
    mpi_comm.cell_parent_mpi_ranks = vector<size_t> { 3, 2, 1, 0, 3, 2 };

    System system;

    // Create a couple of atoms in the first list, owned by rank 3
    CellList list0 {2, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list1 {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list2 {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list3 {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list4 {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list5 {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};

    // Neighbouring cell indices:
    // The atoms in list0 will be sent from rank 3 to ranks 2 and 1
    list0.to_neighbours = vector<size_t> { 1, 2 }; // rank 3 -> 2, 1
    list1.to_neighbours = vector<size_t> { 2, 3 };
    list2.to_neighbours = vector<size_t> { 3, 4 };
    list3.to_neighbours = vector<size_t> { 4, 5 };
    list4.to_neighbours = vector<size_t> { 5, 0 };
    list5.to_neighbours = vector<size_t> { 0, 1 };

    // This vector should match the positions
    const vector<real> xs_expected {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    const vector<real> fs_expected (xs_expected.size(), 0.0);

    // Add the atoms first to the owning rank
    if (is_rank(3, mpi_comm))
    {
        list0.add_atom(0.0, 1.0, 2.0);
        list0.add_atom(3.0, 4.0, 5.0);

        ASSERT_EQ_VEC(list0.xs, xs_expected,
            "the initial positions in the parent rank are incorrect");
    }
    // And to the receivers, for it to be replaced
    else if (is_rank(2, mpi_comm))
    {
        list0.add_atom(0.1, 0.2, 0.3);

        // Old forces, too
        list0.fs = vector<real> {0.4, 0.5, 0.6};
    }
    else if (is_rank(1, mpi_comm))
    {
        list0.add_atom(0.1, 0.2, 0.3);
        list0.add_atom(1.1, 1.2, 1.3);
        list0.add_atom(2.1, 2.2, 2.3);
        list0.add_atom(3.1, 3.2, 3.3);
        list0.add_atom(4.1, 4.2, 4.3);

        list0.fs = vector<real> {
            0.0, 0.2, 0.4,
            0.6, 0.8, 1.0,
            1.2, 1.4, 1.4,
            1.6, 1.8, 2.0,
            2.2, 2.4, 2.6,
        };
    }
    else
    {
        ASSERT_EQ(list0.xs.size(), 0,
            "other ranks do not start with zeroed positions");
    }

    system.cell_lists.push_back(list0);
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);
    system.cell_lists.push_back(list3);
    system.cell_lists.push_back(list4);
    system.cell_lists.push_back(list5);

    // Setup the communication records
    mpi_create_cell_comm_groups(mpi_comm, system.cell_lists);
    mpi_create_sending_and_receiving_cell_lists_for_ranks(
        mpi_comm, system.cell_lists
    );

    // Transfer the atoms and verify that they (the positions) arrived
    mpi_synchronize_interaction_cell_lists(system, mpi_comm);

    // The sender should still have them
    if (is_rank(3, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ_VEC(list.xs, xs_expected,
            "the final positions in the parent rank are incorrect");
    }

    // As should rank 2
    if (is_rank(2, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ_VEC(list.xs, xs_expected,
            "rank 2 did not receive the expected positions");
        ASSERT_EQ_VEC(list.fs, fs_expected,
            "rank 2 did not created nulled-memory for the forces");
    }

    // ... and rank 1
    if (is_rank(1, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ_VEC(list.xs, xs_expected,
            "rank 1 did not receive the expected positions");
        ASSERT_EQ_VEC(list.fs, fs_expected,
            "rank 1 did not created nulled-memory for the forces");
    }

    // The remainining rank 0 should still not have any
    if (is_rank(0, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ(list0.xs.size(), 0,
            "rank 0 received positions, but it shouldn't had");
        ASSERT_EQ(list0.fs.size(), 0,
            "rank 0 received forces, but it shouldn't had");
    }

    // All the remaining cell lists should be empty, for all ranks
    for (auto i = 1; i < system.cell_lists.size(); ++i)
    {
        const auto& list = system.cell_lists.at(i);
        ASSERT_EQ(list.xs.size(), 0,
            "a rank has unexpected positions in a cell that should be empty");
        ASSERT_EQ(list.vs.size(), 0,
            "a rank has unexpected velocities in a cell that should be empty");
        ASSERT_EQ(list.fs.size(), 0,
            "a rank has unexpected forces in a cell that should be empty");
    }
)

ADD_TEST(test_transmitting_empty_cells_erases_cell_data,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // Same setup as the above test
    // 4 ranks with a few cells each
    mpi_comm.cell_parent_mpi_ranks = vector<size_t> { 3, 2, 1, 0, 3, 2 };

    System system;

    // Create a couple of atoms in the first list, owned by rank 3
    CellList list0 {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list1 {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list2 {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list3 {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list4 {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list5 {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};

    // Neighbouring cell indices:
    // The atoms in list0 will be sent from rank 3 to ranks 2 and 1
    list0.to_neighbours = vector<size_t> { 1, 2 }; // rank 3 -> 2, 1
    list1.to_neighbours = vector<size_t> { 2, 3 };
    list2.to_neighbours = vector<size_t> { 3, 4 };
    list3.to_neighbours = vector<size_t> { 4, 5 };
    list4.to_neighbours = vector<size_t> { 5, 0 };
    list5.to_neighbours = vector<size_t> { 0, 1 };

    // This vector should match the positions
    const vector<real> xs_expected {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    const vector<real> fs_expected (xs_expected.size(), 0.0);

    // The receivers have old data which will be erased with an empty cell
    if (is_rank(2, mpi_comm))
    {
        list0.add_atom(0.1, 0.2, 0.3);

        // Old forces, too
        list0.fs = vector<real> {0.4, 0.5, 0.6};
    }
    else if (is_rank(1, mpi_comm))
    {
        list0.add_atom(0.1, 0.2, 0.3);
        list0.add_atom(1.1, 1.2, 1.3);
        list0.add_atom(2.1, 2.2, 2.3);
        list0.add_atom(3.1, 3.2, 3.3);
        list0.add_atom(4.1, 4.2, 4.3);

        list0.fs = vector<real> {
            0.0, 0.2, 0.4,
            0.6, 0.8, 1.0,
            1.2, 1.4, 1.4,
            1.6, 1.8, 2.0,
            2.2, 2.4, 2.6,
        };
    }
    else
    {
        ASSERT_EQ(list0.xs.size(), 0,
            "other ranks do not start with zeroed positions");
    }

    system.cell_lists.push_back(list0);
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);
    system.cell_lists.push_back(list3);
    system.cell_lists.push_back(list4);
    system.cell_lists.push_back(list5);

    // Setup the communication records
    mpi_create_cell_comm_groups(mpi_comm, system.cell_lists);
    mpi_create_sending_and_receiving_cell_lists_for_ranks(
        mpi_comm, system.cell_lists
    );

    // Transfer the atoms and verify that they (the positions) arrived
    mpi_synchronize_interaction_cell_lists(system, mpi_comm);

    if (is_rank(2, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ(list.num_atoms(), 0,
            "rank 2 did not erase its atoms");
        ASSERT_EQ(list.xs.size(), 0,
            "rank 2 did not erase its atoms");
        ASSERT_EQ(list.fs.size(), 0,
            "rank 2 did not erase its atoms");
    }

    // ... and rank 1
    if (is_rank(1, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ(list.num_atoms(), 0,
            "rank 1 did not erase its atoms");
        ASSERT_EQ(list.xs.size(), 0,
            "rank 1 did not erase its atoms");
        ASSERT_EQ(list.fs.size(), 0,
            "rank 1 did not erase its atoms");
    }
)

ADD_TEST(test_collect_forces_from_transmitted_cell_lists,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // Same setup as the transmission tests
    // 4 ranks with a few cells each
    mpi_comm.cell_parent_mpi_ranks = vector<size_t> { 3, 2, 1, 0, 3, 2 };

    System system;

    // Create a couple of atoms in the first list, owned by rank 3
    CellList list0 {2, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list1 {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list2 {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list3 {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list4 {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    CellList list5 {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};

    // Neighbouring cell indices:
    // The atoms in list0 will be sent from rank 3 to ranks 2 and 1
    list0.to_neighbours = vector<size_t> { 1, 2 }; // ranks 2, 1 -> 3
    list1.to_neighbours = vector<size_t> { 2, 3 }; // ranks 1, 0 -> 2
    list2.to_neighbours = vector<size_t> { 3, 4 }; // ranks 0, 3 -> 1
    list3.to_neighbours = vector<size_t> { 4, 5 }; // ranks 3, 2 -> 0
    list4.to_neighbours = vector<size_t> { 5, 0 }; // ranks 2, 3 -> 3
    list5.to_neighbours = vector<size_t> { 0, 1 }; // ranks 3, 2 -> 2

    // All ranks have the same number of atoms this time: only those from
    // 2, 1 will be collected into rank 3 though
    list0.add_atom(0.0, 1.0, 2.0);
    list0.add_atom(3.0, 4.0, 5.0);

    // Some forces to ranks 3, 2 and 1:
    // use the rank number as a multiplier to change them slightly
    if (is_rank(3, mpi_comm) || is_rank(2, mpi_comm) || is_rank(1, mpi_comm))
    {

        const auto rank = mpi_comm.rank;
        list0.fs = vector<real> {
            rank * 0.1, rank * 0.2, rank * 0.3,
            rank * 0.4, rank * 0.5, rank * 0.6
        };

        ASSERT_EQ(list0.fs.size(), list0.xs.size(),
            "setup error: wrong number of forces set for a transmitting rank");
    }
    // Some other forces for the un-synced rank cell0
    else
    {
        list0.fs = vector<real> { 3.0, 5.0, 7.0, 11.0, 13.0, 17.0 };
        ASSERT_EQ(list0.fs.size(), list0.xs.size(),
            "setup error: wrong number of forces set for rank 0");
    }

    // This vector should match the summed forces in cell0:
    // the sum of the rank indices is 6
    const vector<real> fs_expected {
        6 * 0.1, 6 * 0.2, 6 * 0.3, 6 * 0.4, 6 * 0.5, 6 * 0.6
    };

    system.cell_lists.push_back(list0);
    system.cell_lists.push_back(list1);
    system.cell_lists.push_back(list2);
    system.cell_lists.push_back(list3);
    system.cell_lists.push_back(list4);
    system.cell_lists.push_back(list5);

    // Setup the communication records
    mpi_create_cell_comm_groups(mpi_comm, system.cell_lists);
    mpi_create_sending_and_receiving_cell_lists_for_ranks(
        mpi_comm, system.cell_lists
    );

    // Get the unchanged force vectors for the rank
    const vector<real> fs_init(list0.fs.cbegin(), list0.fs.cend());

    // Transfer the atoms and verify that they (the positions) arrived
    mpi_collect_forces_from_interaction_cell_lists(system, mpi_comm);

    // The parent rank should have the sum
    if (is_rank(3, mpi_comm))
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ_VEC(list.fs, fs_expected,
            "the final forces in the owning rank 3 are incorrect");
    }
    // All others should have their initial values
    else
    {
        const auto& list = system.cell_lists[0];
        ASSERT_EQ_VEC(list.fs, fs_init,
            "the final forces in the remaining ranks are incorrect");
    }

    for (auto i = 1; i < system.cell_lists.size(); ++i)
    {
        const auto& list = system.cell_lists[i];
        ASSERT_EQ(list.fs.size(), 0, "a cell got unexpected forces");
    }
)

ADD_TEST(test_reset_received_cells,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // All ranks have atoms in all cells. Only those belonging to each
    // rank will remain after resetting the received lists.
    mpi_comm.mpi_rank_received_cells = vector<vector<size_t>> {
        vector<size_t> { 5, 4 },
        vector<size_t> { 3, 2 },
        vector<size_t> { 1 },
        vector<size_t> { 0 }
    };

    // Create a couple of atoms in the first list, owned by rank 3
    System system;

    system.cell_lists = vector<CellList> {
        CellList {2, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {2, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {2, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {2, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {2, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {2, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}}
    };

    for (auto& list : system.cell_lists)
    {
        list.add_atom(0.0, 1.0, 2.0);
        list.add_atom(3.0, 4.0, 5.0);
    }

    reset_received_cell_lists(system, mpi_comm);

    const auto& cell_lists = system.cell_lists;

    for (auto rank = 0; rank < 4; ++rank)
    {
        for (auto i = 0; i < cell_lists.size(); ++i)
        {
            const auto& list = cell_lists.at(i);

            if (is_rank(rank, mpi_comm))
            {
                if (in_set(
                    static_cast<size_t>(i),
                    mpi_comm.mpi_rank_received_cells.at(rank)
                ))
                {
                    ASSERT_EQ(list.num_atoms(), 0,
                        "a rank was not cleared of atoms in a received cell");
                }
                else
                {
                    ASSERT_EQ(list.num_atoms(), 2,
                        "a rank lost cells in an owned rank");
                }
            }
        }
    }
)

ADD_TEST(test_get_cells_to_transmit_to_other_ranks,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // Construct a system with:
    // 6 cells, divided onto 4 ranks
    vector<CellList> cell_lists {
        CellList {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}}
    };

    mpi_comm.mpi_rank_owned_cells = vector<vector<size_t>> {
        vector<size_t> { 0, 1 },
        vector<size_t> { 2, 3 },
        vector<size_t> { 4 },
        vector<size_t> { 5 }
    };

    mpi_comm.mpi_rank_non_owned_cells = vector<vector<size_t>> {
        vector<size_t> { 2, 3, 4, 5 },
        vector<size_t> { 0, 1, 4, 5 },
        vector<size_t> { 0, 1, 2, 3, 5 },
        vector<size_t> { 0, 1, 2, 3, 4 }
    };

    mpi_comm.cell_parent_mpi_ranks = get_cell_mpi_ranks(
        mpi_comm.mpi_rank_owned_cells
    );

    // Ranks 0 and 2 will transmit atoms to cell 2 and 3 of rank 1
    if (is_rank(0, mpi_comm))
    {
        cell_lists.at(2).add_atom(10.0, 11.0, 12.0);
        cell_lists.at(3).add_atom(10.0, 11.0, 12.0);
    }
    else if (is_rank(2, mpi_comm))
    {
        cell_lists.at(2).add_atom(6.0, 7.0, 8.0);
        cell_lists.at(3).add_atom(6.0, 7.0, 8.0);
    }
    // Rank 3 sends from cell 0 to rank 0
    else if (is_rank(3, mpi_comm))
    {
        cell_lists.at(0).add_atom(6.0, 7.0, 8.0);
    }

    const auto cells_to_transmit_per_rank = get_cells_to_transmit(
        cell_lists, mpi_comm
    );

    for (size_t rank = 0; rank < 4; ++rank)
    {
        if (is_rank(0, mpi_comm) && (rank == 1))
        {
            const auto expected = vector<size_t> { 2, 3 };
            ASSERT_EQ_VEC(cells_to_transmit_per_rank.at(rank), expected,
                "rank 0 does not transmit the correct cells to rank 1");
        }
        else if (is_rank(2, mpi_comm) && (rank == 1))
        {
            const auto expected = vector<size_t> { 2, 3 };
            ASSERT_EQ_VEC(cells_to_transmit_per_rank.at(rank), expected,
                "rank 2 does not transmit the correct cells to rank 1");
        }
        else if (is_rank(3, mpi_comm) && (rank == 0))
        {
            const auto expected = vector<size_t> { 0 };
            ASSERT_EQ_VEC(cells_to_transmit_per_rank.at(rank), expected,
                "rank 3 does not transmit the correct cells to rank 0");
        }
        else
        {
            ASSERT_EQ(cells_to_transmit_per_rank.at(rank).size(), 0,
                "a rank will transmit unexpected cells");
        }
    }
)

ADD_TEST(test_sync_num_of_cells_to_receive_from_other_ranks,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // Construct a system with:
    // 6 cells, divided onto 4 ranks
    vector<CellList> cell_lists {
        CellList {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}}
    };

    mpi_comm.mpi_rank_owned_cells = vector<vector<size_t>> {
        vector<size_t> { 0, 1 },
        vector<size_t> { 2, 3 },
        vector<size_t> { 4 },
        vector<size_t> { 5 }
    };

    mpi_comm.mpi_rank_non_owned_cells = vector<vector<size_t>> {
        vector<size_t> { 2, 3, 4, 5 },
        vector<size_t> { 0, 1, 4, 5 },
        vector<size_t> { 0, 1, 2, 3, 5 },
        vector<size_t> { 0, 1, 2, 3, 4 }
    };

    mpi_comm.cell_parent_mpi_ranks = get_cell_mpi_ranks(
        mpi_comm.mpi_rank_owned_cells
    );

    // Ranks 0 and 2 will transmit atoms to cell 2 and 3 of rank 1
    if (is_rank(0, mpi_comm))
    {
        cell_lists.at(2).add_atom(10.0, 11.0, 12.0);
        cell_lists.at(3).add_atom(10.0, 11.0, 12.0);
    }
    else if (is_rank(2, mpi_comm))
    {
        cell_lists.at(2).add_atom(6.0, 7.0, 8.0);
        cell_lists.at(3).add_atom(6.0, 7.0, 8.0);
    }
    // Rank 3 sends from cell 0 to rank 0 and cell 4 to 2
    else if (is_rank(3, mpi_comm))
    {
        cell_lists.at(0).add_atom(6.0, 7.0, 8.0);
        cell_lists.at(4).add_atom(6.0, 7.0, 8.0);
    }

    const auto cells_to_transmit_per_rank = get_cells_to_transmit(
        cell_lists, mpi_comm
    );
    const auto num_cells_to_receive_per_rank
        = mpi_sync_number_of_transmitted_cells(
            cells_to_transmit_per_rank, mpi_comm);

    if (is_rank(0, mpi_comm))
    {
        const vector<uint64_t> expected {0, 0, 0, 1};
        ASSERT_EQ_VEC(num_cells_to_receive_per_rank, expected,
            "rank 0 does not receive the correct number of cells"
        );
    }
    else if (is_rank(1, mpi_comm))
    {
        const vector<uint64_t> expected {2, 0, 2, 0};
        ASSERT_EQ_VEC(num_cells_to_receive_per_rank, expected,
            "rank 0 does not receive the correct number of cells"
        );
    }
    else if (is_rank(2, mpi_comm))
    {
        const vector<uint64_t> expected {0, 0, 0, 1};
        ASSERT_EQ_VEC(num_cells_to_receive_per_rank, expected,
            "rank 0 does not receive the correct number of cells"
        );
    }
    else if (is_rank(3, mpi_comm))
    {
        const vector<uint64_t> expected {0, 0, 0, 0};
        ASSERT_EQ_VEC(num_cells_to_receive_per_rank, expected,
            "rank 0 does not receive the correct number of cells"
        );
    }

)

ADD_TEST(test_move_atoms_to_new_owning_ranks,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // Construct a system with:
    // 6 cells, divided onto 4 ranks
    vector<CellList> cell_lists {
        CellList {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}}
    };

    mpi_comm.mpi_rank_owned_cells = vector<vector<size_t>> {
        vector<size_t> { 0, 1 },
        vector<size_t> { 2, 3 },
        vector<size_t> { 4 },
        vector<size_t> { 5 }
    };

    mpi_comm.mpi_rank_non_owned_cells = vector<vector<size_t>> {
        vector<size_t> { 2, 3, 4, 5 },
        vector<size_t> { 0, 1, 4, 5 },
        vector<size_t> { 0, 1, 2, 3, 5 },
        vector<size_t> { 0, 1, 2, 3, 4 }
    };

    mpi_comm.cell_parent_mpi_ranks = get_cell_mpi_ranks(
        mpi_comm.mpi_rank_owned_cells
    );

    // Rank 0 and rank 2 have atoms in cell 2, which belongs to rank 1.
    if (is_rank(0, mpi_comm))
    {
        cell_lists.at(2).add_atom(0.0, 1.0, 2.0);
        cell_lists.at(2).add_atom(3.0, 4.0, 5.0);
        cell_lists.at(2).vs = vector<real> {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};

        // Forces will not be transmitted
        cell_lists.at(2).fs = vector<real> {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
    }
    else if (is_rank(2, mpi_comm))
    {
        cell_lists.at(2).add_atom(6.0, 7.0, 8.0);
        cell_lists.at(2).vs = vector<real> {0.6, 0.7, 0.8};
        cell_lists.at(2).fs = vector<real> {0.6, 0.7, 0.8};
    }
    // They should be merged onto rank 1's list, which already has an atom
    else if (is_rank(1, mpi_comm))
    {
        cell_lists.at(2).add_atom(10.0, 11.0, 12.0);
        cell_lists.at(2).vs = vector<real> {10.1, 11.1, 12.1};

        // Only these forces will be kept
        cell_lists.at(2).fs = vector<real> {10.1, 11.1, 12.1};
    }

    System system;
    system.cell_lists = move(cell_lists);

    mpi_move_atoms_to_owning_ranks(system, mpi_comm);

    // Check rank 1's atoms
    if (is_rank(1, mpi_comm))
    {
        for (size_t i = 0; i < cell_lists.size(); ++i)
        {
            const auto& list = system.cell_lists.at(i);

            if (i == 2)
            {
                ASSERT_EQ(list.num_atoms(), 4,
                    "rank 1 does not have the correct number of atoms "
                    "in cell 2");

                const vector<real> xs_expected {
                        10.0, 11.0, 12.0,
                        0.0, 1.0, 2.0,
                        3.0, 4.0, 5.0,
                        6.0, 7.0, 8.0
                };

                const vector<real> vs_expected {
                    10.1, 11.1, 12.1,
                    0.0, 0.1, 0.2,
                    0.3, 0.4, 0.5,
                    0.6, 0.7, 0.8
                };

                const vector<real> fs_expected {
                    10.1, 11.1, 12.1,
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0
                };

                ASSERT_EQ_VEC(list.xs, xs_expected,
                    "rank 1 does not have the correct positions in cell 1 ");
                ASSERT_EQ_VEC(list.vs, vs_expected,
                    "rank 1 does not have the correct velocities in cell 1 ");
                ASSERT_EQ_VEC(list.fs, fs_expected,
                    "rank 1 does not have the correct forces in cell 1 ");
            }
            else
            {
                ASSERT_EQ(list.num_atoms(), 0,
                    "rank 1 has atoms in other cells");
            }
        }
    }
)

ADD_TEST(test_collect_all_data_on_master_rank,
    MPIRank mpi_comm;
    init_MPI(mpi_comm);

    // Construct a system with:
    // 6 cells, divided onto 4 ranks
    System system;

    vector<CellList> cell_lists {
        CellList {0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 0.0, 0.0}, {1.0, 1.0, 1.0}},
        CellList {0, {2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}}
    };

    mpi_comm.mpi_rank_owned_cells = vector<vector<size_t>> {
        vector<size_t> { 0, 1 },
        vector<size_t> { 2, 3 },
        vector<size_t> { 4 },
        vector<size_t> { 5 }
    };

    mpi_comm.mpi_rank_non_owned_cells = get_non_owned_cells_per_rank(
        mpi_comm.mpi_rank_owned_cells, 6
    );
    mpi_comm.cell_parent_mpi_ranks = get_cell_mpi_ranks(
        mpi_comm.mpi_rank_owned_cells
    );

    // All ranks have two atoms in their cells
    if (is_rank(0, mpi_comm))
    {
        cell_lists.at(0).add_atom(0.1, 0.2, 0.3);
        cell_lists.at(0).add_atom(0.4, 0.5, 0.6);
        cell_lists.at(0).vs = vector<real> {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};

        cell_lists.at(1).add_atom(1.1, 1.2, 1.3);
        cell_lists.at(1).add_atom(1.4, 1.5, 1.6);
        cell_lists.at(1).vs = vector<real> {1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
    }
    else if (is_rank(1, mpi_comm))
    {
        cell_lists.at(2).add_atom(2.1, 2.2, 2.3);
        cell_lists.at(2).add_atom(2.4, 2.5, 2.6);
        cell_lists.at(2).vs = vector<real> {2.0, 2.1, 2.2, 2.3, 2.4, 2.5};

        cell_lists.at(3).add_atom(3.1, 3.2, 3.3);
        cell_lists.at(3).add_atom(3.4, 3.5, 3.6);
        cell_lists.at(3).vs = vector<real> {3.0, 3.1, 3.2, 3.3, 3.4, 3.5};
    }
    else if (is_rank(2, mpi_comm))
    {
        cell_lists.at(4).add_atom(4.1, 4.2, 4.3);
        cell_lists.at(4).add_atom(4.4, 4.5, 4.6);
        cell_lists.at(4).vs = vector<real> {4.0, 4.1, 4.2, 4.3, 4.4, 4.5};
    }
    else if (is_rank(3, mpi_comm))
    {
        cell_lists.at(5).add_atom(5.1, 5.2, 5.3);
        cell_lists.at(5).add_atom(5.4, 5.5, 5.6);
        cell_lists.at(5).vs = vector<real> {5.0, 5.1, 5.2, 5.3, 5.4, 5.5};
    }

    system.cell_lists = move(cell_lists);
    mpi_collect_atoms_to_master(system, mpi_comm);

    if (is_master(mpi_comm))
    {
        ASSERT_EQ(system.num_atoms(), 12,
            "master rank does not have all atoms after sync");

        const auto xs_cell0 = vector<real> {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
        const auto vs_cell0 = vector<real> {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};

        const auto xs_cell1 = vector<real> {1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
        const auto vs_cell1 = vector<real> {1.0, 1.1, 1.2, 1.3, 1.4, 1.5};

        const auto xs_cell2 = vector<real> {2.1, 2.2, 2.3, 2.4, 2.5, 2.6};
        const auto vs_cell2 = vector<real> {2.0, 2.1, 2.2, 2.3, 2.4, 2.5};

        const auto xs_cell3 = vector<real> {3.1, 3.2, 3.3, 3.4, 3.5, 3.6};
        const auto vs_cell3 = vector<real> {3.0, 3.1, 3.2, 3.3, 3.4, 3.5};

        const auto xs_cell4 = vector<real> {4.1, 4.2, 4.3, 4.4, 4.5, 4.6};
        const auto vs_cell4 = vector<real> {4.0, 4.1, 4.2, 4.3, 4.4, 4.5};

        const auto xs_cell5 = vector<real> {5.1, 5.2, 5.3, 5.4, 5.5, 5.6};
        const auto vs_cell5 = vector<real> {5.0, 5.1, 5.2, 5.3, 5.4, 5.5};

        const auto lists = system.cell_lists;

        ASSERT_EQ_VEC(lists[0].xs, xs_cell0,
            "cell 0 does not have the correct positions");
        ASSERT_EQ_VEC(lists[0].vs, vs_cell0,
            "cell 0 does not have the correct velocities");
        ASSERT_EQ_VEC(lists[1].xs, xs_cell1,
            "cell 1 does not have the correct positions");
        ASSERT_EQ_VEC(lists[1].vs, vs_cell1,
            "cell 1 does not have the correct velocities");
        ASSERT_EQ_VEC(lists[2].xs, xs_cell2,
            "cell 2 does not have the correct positions");
        ASSERT_EQ_VEC(lists[2].vs, vs_cell2,
            "cell 2 does not have the correct velocities");
        ASSERT_EQ_VEC(lists[3].xs, xs_cell3,
            "cell 3 does not have the correct positions");
        ASSERT_EQ_VEC(lists[3].vs, vs_cell3,
            "cell 3 does not have the correct velocities");
        ASSERT_EQ_VEC(lists[4].xs, xs_cell4,
            "cell 4 does not have the correct positions");
        ASSERT_EQ_VEC(lists[4].vs, vs_cell4,
            "cell 4 does not have the correct velocities");
        ASSERT_EQ_VEC(lists[5].xs, xs_cell5,
            "cell 5 does not have the correct positions");
        ASSERT_EQ_VEC(lists[5].vs, vs_cell5,
            "cell 5 does not have the correct velocities");
    }
)

RUN_TESTS(
    MPI_Init(0, nullptr);
    test_init_mpi_for_test_suite();

    // Rank-cell-ownership data setup
    test_cell_indices_are_correctly_divided();
    test_get_cell_ranks_from_owned_cells_per_rank_lists();
    test_get_owned_cells_per_rank();
    test_get_owned_cells_per_rank_even_dividers();
    test_get_owned_cells_per_rank_odd_dividers();
    test_get_non_owned_cells_for_ranks();
    test_get_cell_ownership_per_rank_metadata();
    test_divide_cell_lists_onto_proper_ranks();

    // Creation of MPI required data and records
    test_create_all_cells_for_all_mpi_ranks_with_correct_fields();
    test_create_cell_mpi_communication_groups_for_sending();
    test_create_record_of_which_cells_every_rank_will_receive();
    test_create_record_of_which_cells_every_rank_will_send();
    test_sending_and_receiving_cell_lists_are_consistent();
    test_fill_in_communicators_and_sendrecv_cell_lists();

    // Syncing of data between threads
    test_sync_num_of_transmitted_positions();
    test_reserve_memory_for_recv_cell_lists_positions_and_forces();
    test_transmitting_a_few_atom_positions_works();
    test_transmitting_positions_erases_old_data_on_receivers();
    test_transmitting_empty_cells_erases_cell_data();
    test_collect_forces_from_transmitted_cell_lists();
    test_get_cells_to_transmit_to_other_ranks();
    test_sync_num_of_cells_to_receive_from_other_ranks();
    test_move_atoms_to_new_owning_ranks();
    test_collect_all_data_on_master_rank();

    // Misc
    test_reset_received_cells();

    MPI_Finalize();
)
