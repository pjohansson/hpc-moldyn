#include <ctime>
#include <iomanip>
#include <iostream>

#include "analytics.h"
#include "conf.h"
#include "integrator.h"
#include "params.h"

struct InputArgs {
    std::string input_conf,
                output_conf;
    uint64_t num_steps;
};

enum InputArgPosition {
    CONF = 1,
    OUTPUT = 2,
    NSTEPS = 3,
    NARGS = 4
};

static bool read_cli_arguments(InputArgs& input_args,
                               const int argc,
                               const char* argv[])
{
    if (argc < NARGS)
    {
        return false;
    }

    input_args.input_conf = static_cast<std::string>(argv[CONF]);
    input_args.output_conf = static_cast<std::string>(argv[OUTPUT]);
    input_args.num_steps = static_cast<uint64_t>(
        std::strtoull(argv[NSTEPS], nullptr, 10)
    );

    return true;
}

int main(const int argc, const char* argv[])
{
    InputArgs input_args;
    if (!read_cli_arguments(input_args, argc, argv))
    {
        std::cerr << "usage: " << argv[0] << " CONF OUTPUT NSTEPS\n"
                  << "\n"
                  << "where\n"
                  << "  CONF is an input configuration in .gro format\n"
                  << "  OUTPUT is the final output configuration\n"
                  << "  NSTEPS is the number of steps to run\n";
        return 1;
    }

    std::cout << "Input arguments: "
              << "conf = " << input_args.input_conf
              << ", output = " << input_args.output_conf
              << ", num_steps = " << input_args.num_steps
              << "\n\n";

    std::cerr << "Reading configuration ... ";
    auto system = read_conf_from_grofile(input_args.input_conf, DefaultFF.sigma);
    std::cerr << "done.\n";

    std::cerr << "Setting up system ... ";
    create_cell_lists(system, DefaultFF.rcut);
    std::cerr << "done.\n";

    std::cerr << '\n';
    describe_system_config(system);
    std::cerr << '\n';

    Benchmark benchmark;
    Energetics energy;

    std::cerr << "Simulating " << input_args.num_steps << " steps:\n";

    for (unsigned step = 0; step < input_args.num_steps; ++step)
    {
        if ((step % 10) == 0)
        {
            std::cerr << "\rstep " << step;
        }

        run_velocity_verlet(system, benchmark, DefaultFF, DefaultOpts);

        if ((step % DefaultOpts.energy_calc) == 0)
        {
            calculate_system_energetics(energy, system, DefaultFF);
        }
    }

    benchmark.finalize();

    std::cerr << "\r                                                  \r"
              << "Finished.\n\n";

    std::cerr << "Writing final system to disk as '"
        << input_args.output_conf
        << "' ... ";
    write_conf_to_grofile(system, input_args.output_conf, DefaultFF.sigma);
    std::cerr << "done.\n";

    std::cerr << '\n';
    print_energetics(energy);
    std::cerr << '\n';
    print_benchmark(benchmark);

    auto t = std::time(nullptr);
    std::cerr << "\nFinished simulation at "
        << std::put_time(std::localtime(&t), "%c") << ".\n";

    return 0;
}
