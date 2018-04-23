#include <cctype>
#include <exception>
#include <fstream>
#include <iostream>
#include <utility>

#include "params.h"

template <typename Iterator>
static Iterator find_string_range(Iterator& it, const Iterator it_end)
{
        while (isblank(*it) && it != it_end)
        {
            ++it;
        }

        const auto string_begin = it++;

        while ((isalnum(*it) || (*it == '.')) && it != it_end)
        {
            ++it;
        }

        return string_begin;
}

bool read_parameter_file(const std::string path, Options& opts)
{
    std::ifstream ifs { path, std::ifstream::in };

    if (!ifs)
    {
        std::cerr << "error: could not open parameter file at '" << path
            << "' for writing" << '\n';

        return false;
    }

    constexpr size_t buflen = 256;
    std::string buffer (buflen, ' ');

    while (ifs.good())
    {
        getline(ifs, buffer);

        auto it = buffer.cbegin();
        const auto field_begin = find_string_range(it, buffer.cend());

        if (*field_begin == '#')
        {
            continue;
        }

        auto field = std::string(field_begin, it);
        for (auto& c : field)
        {
            c = std::toupper(c);
        }

        const auto value_begin = find_string_range(it, buffer.cend());
        auto value = std::string(value_begin, it);
        for (auto& c : value)
        {
            c = std::toupper(c);
        }

        try {
            if (field == "NSTEPS")
            {
                opts.num_steps = static_cast<unsigned>(std::stoull(value));
            }
            else if (field == "DT")
            {
                const auto dt = static_cast<double>(std::stod(value));
                opts.set_dt(dt);
            }
            else if (field == "ENERGYOUT")
            {
                opts.energy_calc = static_cast<unsigned>(std::stoull(value));
            }
            else if (field == "TRAJOUT")
            {
                opts.traj_stride = static_cast<unsigned>(std::stoull(value));
            }
            else if (field == "GENVEL")
            {
                if (value == "TRUE")
                {
                    opts.gen_velocities = true;
                }
                else if (value != "FALSE")
                {
                    std::cerr
                        << "error: "
                        << "'GENVEL' must be either 'TRUE/FALSE', was '"
                        << value << '\'' << '\n';

                    return false;
                }
            }
            else if (field == "GENTEMP")
            {
                opts.gen_temp = static_cast<double>(std::stod(value));
            }
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr
                << "error: '" << field
                << "' got a bad value '" << value << '\'' << '\n';

            return false;
        }
    }

    return true;
}
