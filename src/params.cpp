#include <cctype>
#include <exception>
#include <fstream>
#include <iostream>
#include <utility>

#include "params.h"

// Advance the input iterator to the end of the input string, indicated
// by the first non-leading non-alphanumeric character (including '+', '-', '.')
// Return the iterator position of the first non-whitespace character to trim
// the string.
template <typename Iterator>
static Iterator find_string_range(Iterator& it, const Iterator it_end)
{
        while (isblank(*it) && it != it_end)
        {
            ++it;
        }

        const auto string_begin = it++;

        while (
            (isalnum(*it) || (*it == '.') || (*it == '-') || (*it == '+'))
            && (it != it_end)
        )
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

ForceField read_forcefield_file(const std::string path)
{
    std::ifstream ifs { path, std::ifstream::in };

    if (!ifs)
    {
        std::cerr << "error: could not open force field file at '" << path
            << "' for writing" << '\n';

        return InvalidFF;
    }

    constexpr size_t buflen = 256;
    std::string buffer (buflen, ' ');

    // Uninitialized values which will not be used, since we check that
    // all values are read from the file before returning.
    real epsilon       = -1.0,
         sigma         = -1.0,
         mass          = -1.0,
         rcut          = -1.0,
         wall_constant = -1.0;
    bool bEpsilon = false,
         bSigma = false,
         bMass = false,
         bRcut = false,
         bWallConstant = false;

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
        const auto value = std::string(value_begin, it);

        try {
            if (field == "EPSILON")
            {
                epsilon = static_cast<real>(std::stod(value));
                bEpsilon = true;
            }
            else if (field == "SIGMA")
            {
                sigma = static_cast<real>(std::stod(value));
                bSigma = true;
            }
            else if (field == "MASS")
            {
                mass = static_cast<real>(std::stod(value));
                bMass = true;
            }
            else if (field == "CUTOFF")
            {
                rcut = static_cast<real>(std::stod(value));
                bRcut = true;
            }
            else if (field == "WALLCONSTANT")
            {
                wall_constant = static_cast<real>(std::stod(value));
                bWallConstant = true;
            }
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr
                << "error: '" << field
                << "' got a bad value '" << value << '\'' << '\n';

            return InvalidFF;
        }
    }

    // All parameters have to be set.
    if (!(bEpsilon && bSigma && bMass && bRcut && bWallConstant))
    {
        std::string missing_fields;

        if (!bEpsilon)
        {
            missing_fields.append(" EPSILON");
        }
        if (!bSigma)
        {
            missing_fields.append(" SIGMA");
        }
        if (!bMass)
        {
            missing_fields.append(" MASS");
        }
        if (!bRcut)
        {
            missing_fields.append(" CUTOFF");
        }
        if (!bWallConstant)
        {
            missing_fields.append(" WALLCONSTANT");
        }

        std::cerr
            << "error: force field file at '" << path << "' missing fields:"
            << missing_fields << '\n';

        return InvalidFF;
    }

    return ForceField(epsilon, sigma, rcut, mass, wall_constant);
}
