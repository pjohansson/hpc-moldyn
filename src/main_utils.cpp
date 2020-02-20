enum class ParseResult {
    Error,
    Advance,
    Ok,
};

struct InputArgs {
    InputArgs()
    :verbose { false },
     read_forcefield { false },
     param_file { "params.dat" },
     forcefield_file { "" },
     traj_file { "traj.gro" }
    {}

    bool verbose, read_forcefield;
    std::string input_conf,
                output_conf,
                param_file,
                forcefield_file,
                traj_file;
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
    else if (argument == "f")
    {
        if (!parse_string_from_next(it, it_end, input_args.forcefield_file))
        {
            return ParseResult::Error;
        }

        return ParseResult::Advance;
    }
    else if (argument == "t")
    {
        if (!parse_string_from_next(it, it_end, input_args.traj_file))
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

    if (input_args.forcefield_file.size() > 0)
    {
        input_args.read_forcefield = true;
    }

    return true;
}

static bool do_step(const size_t step, const size_t stride)
{
    return (stride > 0) && (step > 0) && (step % stride == 0);
}
