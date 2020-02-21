#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>

#include <mpi.h>

#include "analytics.h"
#include "conf.h"
#include "integrator.h"
#include "main_utils.cpp"
#include "mpi_impl.h"
#include "params.h"

// #define DEBUG 1;
// #define DEBUG_HOLD 1;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    MPIRank mpi_comm;
    if (!init_MPI(mpi_comm))
    {
        std::cerr << "error: could not initialize MPI\n";
        return 1;
    };

    Benchmark benchmark;
    Energetics energy;
    InputArgs input_args;
    Options opts;
    System system;

    ForceField ff = UninitializedFF;

    if (is_master(mpi_comm))
    {
        const std::vector<std::string> str_args (argv + 1, argv + argc);

        if (!read_cli_arguments(str_args, input_args))
        {
            std::cerr
                << "usage: " << argv[0] << "[OPTIONS] <CONF> <OUTPUT>\n"
                << "\n"
                << "  CONF is an input configuration in .gro format\n"
                << "  OUTPUT is the final output configuration\n"
                << "\n"
                << "The following optional flags are available:\n"
                << "  -p <params.dat>  File with the simulation parameters\n"
                << "  -f <[.dat]>      File with force field parameters (default: Argon)\n"
                << "  -t <traj.gro>    Write trajectory to this file\n"
                << "  -v               Be loud and noisy\n";

            return 1;
        }

        std::cout << "Input arguments:" << '\n'
                  << "  conf = " << input_args.input_conf << '\n'
                  << "  output = " << input_args.output_conf << '\n'
                  << "  parameter file = " << input_args.param_file << '\n';
        if (input_args.read_forcefield)
        {
            std::cout
                << "  forcefield file = " << input_args.forcefield_file << '\n';
        }
        std::cout << "\n";

        std::cerr << "Reading parameter file ... ";
        if (!read_parameter_file(input_args.param_file, opts))
        {
            return 1;
        }
        std::cerr << "done.\n";

        if (input_args.read_forcefield)
        {
            std::cerr << "Reading force field file ... ";
            ff = read_forcefield_file(input_args.forcefield_file);
            if (!ff.is_valid)
            {
                return 1;
            }
            std::cerr << "done.\n";
        }
        else
        {
            ff = ArgonFF;
            std::cerr << "No force field was specified. Defaulting to Argon parameters.\n";
        }

        std::cerr << "Reading configuration ... ";
        system = read_conf_from_grofile(
            input_args.input_conf, ff.sigma
        );
        std::cerr << "done.\n";

        std::cerr << "Setting up system ... ";
        create_cell_lists(system, ff.rcut);
        std::cerr << "done.\n\n";

        if (opts.gen_velocities)
        {
            std::cerr << "Generating atom velocities around T = " << opts.gen_temp << "  ... ";
            gen_system_velocities(system, opts.gen_temp);
            std::cerr << "done.\n\n";
        }

        if (input_args.verbose)
        {
            describe_system_config(system, ff);
            std::cerr << '\n';
        }

        // Reopen to clear the trajectory before writing frames into it.
        benchmark.start_traj_output_update();
        std::ofstream reset {
            input_args.traj_file, std::ios::out | std::ios::trunc
        };
        benchmark.stop_traj_output_update();
    }

    benchmark.start_mpi_setup_update();

    if (!sync_options(opts, mpi_comm))
    {
        if (is_master(mpi_comm))
        {
            std::cerr << "error: could not sync data between MPI ranks\n";
        }

        MPI_Finalize();

        return 1;
    };

    if (!sync_forcefield(ff, mpi_comm))
    {
        if (is_master(mpi_comm))
        {
            std::cerr << "error: could not sync force field between MPI ranks\n";
        }

        MPI_Finalize();

        return 1;
    };

    mpi_init_cell_lists_and_transfer_atoms(system, mpi_comm);
    mpi_move_atoms_to_owning_ranks(system, mpi_comm);
    mpi_fill_communication_data(mpi_comm, system);

    benchmark.stop_mpi_setup_update();

    // size_t stepout_stride = 10;
    size_t stepout_stride = 10;

    if (is_master(mpi_comm))
    {
        describe_system_config(system, ff);

        std::cerr << "Force field parameters:\n"
            << "  epsilon       = " << ff.epsilon << " [J]\n"
            << "  sigma         = " << ff.sigma << " [nm]\n"
            << "  mass          = " << ff.mass << " [u]\n"
            << "  cutoff        = " << ff.rcut << " [sigma]\n"
            << "  wall constant = " << ff.wall_constant << " [epsilon]\n\n";

        std::cerr << "Simulating " << opts.num_steps << " steps:\n";
    }

    benchmark.start_simulation_time_update();

    for (uint64_t step = 0; step < opts.num_steps; ++step)
    {
        if (is_master(mpi_comm) && ((step % stepout_stride) == 0))
        {
            std::cerr << "\rstep " << step;

            if (step >= 200)
            {
                stepout_stride = 100;
            }
        }

        run_velocity_verlet(system, benchmark, mpi_comm, ff, opts);

        // Do we need to calculate energetics or output the trajectory
        // this step? If so, collect all atoms on the master rank.
        const bool do_collect_on_master
            = do_step(step, opts.energy_calc)
                || do_step(step, opts.traj_stride);

        if (do_collect_on_master)
        {
            mpi_collect_atoms_to_master(system, mpi_comm);
        }

        benchmark.start_energy_calc_update();
        if (is_master(mpi_comm) && do_step(step, opts.energy_calc))
        {
            calculate_system_energetics(energy, system, ff);
        }
        benchmark.stop_energy_calc_update();

        benchmark.start_traj_output_update();
        if (is_master(mpi_comm) && do_step(step, opts.traj_stride))
        {
            write_conf_to_grofile(
                system, input_args.traj_file, ff.sigma, OutputMode::Append
            );
        }
        benchmark.stop_traj_output_update();

        // Remove atoms from other ranks on master
        if (is_master(mpi_comm))
        {
            for (const auto i : mpi_comm.mpi_rank_non_owned_cells.at(0))
            {
                system.cell_lists.at(i).resize_atom_list(0);
            }
        }

        // TODO: time idle
        MPI_Barrier(MPI_COMM_WORLD);
    }

    benchmark.stop_simulation_time_update();

    mpi_collect_atoms_to_master(system, mpi_comm);

    if (is_master(mpi_comm))
    {
        benchmark.finalize();

        std::cerr << "\r                                                  \r"
                  << "Finished.\n\n";


        std::cerr << "Writing final system to disk as '"
            << input_args.output_conf
            << "' ... ";
        write_conf_to_grofile(system, input_args.output_conf, ff.sigma, OutputMode::Replace);
        std::cerr << "done.\n";

        std::cerr << '\n';
        print_energetics(energy, ff);
        std::cerr << '\n';
        print_benchmark(benchmark, opts.num_steps);

        auto t = std::time(nullptr);
        std::cerr << "\nFinished simulation at "
            << std::put_time(std::localtime(&t), "%c") << ".\n";
    }

    MPI_Finalize();

    return 0;
}
