#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>

#include "analytics.h"
#include "conf.h"
#include "integrator.h"
#include "params.h"

enum class ParseResult {
    Error,
    Advance,
    Ok,
};

struct InputArgs {
    InputArgs()
    :verbose { false },
     param_file { "params.dat" }
    {}

    bool verbose;
    std::string input_conf,
                output_conf,
                param_file;
};

template <typename Iterator>
static bool parse_string_from_next(const Iterator it,
                                   const Iterator it_end,
                                   std::string& ref)
{
    const auto it_value = std::next(it);
    if (it_value == it_end)
    {
        std::cerr << "error: option '" << *it
            << "' requires a value but none was given" << '\n';

        return false;
    }

    ref = *it_value;

    return true;
}

template <typename Iterator>
static ParseResult parse_argument(const Iterator it,
                                  const Iterator it_end,
                                  InputArgs& input_args)
{
    const auto argument = (*it).substr(1);

    if (argument == "p")
    {
        if (!parse_string_from_next(it, it_end, input_args.param_file))
        {
            return ParseResult::Error;
        }

        return ParseResult::Advance;
    }
    else if (argument == "v")
    {
        input_args.verbose = true;
    }
    else
    {
        std::cerr << "error: invalid argument '" << *it << '\'' << '\n';

        return ParseResult::Error;
    }

    return ParseResult::Ok;
}

static bool read_cli_arguments(const std::vector<std::string>& args,
                               InputArgs& input_args)
{
    size_t i_positional = 0;

    auto arg = args.cbegin();
    const auto end = args.cend();

    while (arg != end)
    {
        if ((*arg).front() == '-')
        {
            switch(parse_argument(arg, end, input_args))
            {
                case ParseResult::Error:
                    return false;
                case ParseResult::Advance:
                    ++arg;
                default:
                    break;
            }
        }
        else
        {
            switch (i_positional)
            {
                case 0:
                    input_args.input_conf = *arg;
                    break;

                case 1:
                    input_args.output_conf = *arg;
                    break;

                default:
                    std::cerr << "error: too many arguments specified" << '\n';
                    return false;
            }

            ++i_positional;
        }

        ++arg;
    }

    if (i_positional < 2)
    {
        return false;
    }

    return true;
}

static bool do_step(const size_t step, const size_t stride)
{
    return (stride > 0) && (step > 0) && (step % stride == 0);
}

int main(const int argc, const char* argv[])
{
    const std::vector<std::string> str_args (argv + 1, argv + argc);

    InputArgs input_args;
    if (!read_cli_arguments(str_args, input_args))
    {
        std::cerr << "usage: " << argv[0] << "[OPTIONS] <CONF> <OUTPUT>\n"
                  << "\n"
                  << "  CONF is an input configuration in .gro format\n"
                  << "  OUTPUT is the final output configuration\n"
                  << "\n"
                  << "The following optional flags are available:\n"
                  << "  -p <params.dat>  File with the simulation parameters\n"
                  << "  -v               Be loud and noisy\n";

        return 1;
    }

    std::cout << "Input arguments: "
              << "conf = " << input_args.input_conf
              << ", output = " << input_args.output_conf
              << ", parameter file = " << input_args.param_file
              << "\n\n";

    const auto& ff = ArgonFF;

    std::cerr << "Reading configuration ... ";
    auto system = read_conf_from_grofile(input_args.input_conf, ff.sigma);
    std::cerr << "done.\n";

    std::cerr << "Reading parameter file ... ";
    Options opts;
    if (!read_parameter_file(input_args.param_file, opts))
    {
        return 1;
    }
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

    Benchmark benchmark;
    Energetics energy;

    std::cerr << "Simulating " << opts.num_steps << " steps:\n";

    // Reopen to clear the trajectory before writing frames into it.
    benchmark.start_traj_output_update();
    std::string fntraj { "traj.gro" };
    {
        std::ofstream reset { fntraj, std::ios::out | std::ios::trunc };
    }
    benchmark.stop_traj_output_update();

    size_t stepout_stride = 10;

    for (unsigned step = 0; step < opts.num_steps; ++step)
    {
        if ((step % stepout_stride) == 0)
        {
            std::cerr << "\rstep " << step;

            if (step >= 200)
            {
                stepout_stride = 100;
            }
        }

        run_velocity_verlet(system, benchmark, ff, opts);

        benchmark.start_energy_calc_update();
        if (do_step(step, opts.energy_calc))
        {
            calculate_system_energetics(energy, system, ff);
        }
        benchmark.stop_energy_calc_update();

        benchmark.start_traj_output_update();
        if (do_step(step, opts.traj_stride))
        {
            write_conf_to_grofile(system, fntraj, ff.sigma, OutputMode::Append);
        }
        benchmark.stop_traj_output_update();
    }

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
    print_benchmark(benchmark);

    auto t = std::time(nullptr);
    std::cerr << "\nFinished simulation at "
        << std::put_time(std::localtime(&t), "%c") << ".\n";

    return 0;
}
