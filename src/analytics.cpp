#include "conf.h"
#include <iomanip>

void show_atom_cell_list_distribution(const System& system)
{
    constexpr auto max_hist_length = 25;

    uint64_t max_atoms = 0;
    for (const auto& list : system.cell_lists)
    {
        if (list.num_atoms() > max_atoms)
        {
            max_atoms = list.num_atoms();
        }
    }

    unsigned n = 0;
    std::string stars(max_hist_length, ' ');

    std::cerr << "Atom distribution in cell lists:\n";

    for (const auto& list : system.cell_lists)
    {
        const auto frac = list.num_atoms() / static_cast<float>(max_atoms);
        const auto num_chars = static_cast<unsigned>(frac * max_hist_length);

        stars.replace(stars.cbegin(), stars.cend(), max_hist_length, ' ');
        stars.replace(0, num_chars, num_chars, '*');

        std::cerr << std::setw(4) << n << ": "
                  << stars << " (" << list.num_atoms() << ")\n";
        ++n;
    }
}

void describe_system_config(const System& system)
{
    std::cerr << "System title: " << system.title << '\n'
              << "Size:  "
                << system.box_size[XX] << " x "
                << system.box_size[YY] << " x "
                << system.box_size[ZZ] << " (nm^3)\n"
              << "Number of atoms:  " << system.num_atoms() << '\n'
              << "Cell list size:  "
                << system.cell_size[XX] << " x "
                << system.cell_size[YY] << " x "
                << system.cell_size[ZZ] << " (nm^3)\n"
              << "Cell list configuration:  "
                << system.shape[XX] << " x "
                << system.shape[YY] << " x "
                << system.shape[ZZ] << '\n'
              << '\n';

    show_atom_cell_list_distribution(system);
}
