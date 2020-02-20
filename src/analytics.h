#include <chrono>

#include "conf.h"
#include "params.h"

#ifndef ANALYTICS_H
#define ANALYTICS_H

using ChronoTime = decltype(std::chrono::system_clock::now());

class Benchmark {
public:
    // Initialize all timers to 0 and start the total simulation timer.
    Benchmark(void)
    :cell_list_update { 0.0 },
     force_update { 0.0 },
     force_wall_update { 0.0 },
     position_update { 0.0 },
     velocity_update { 0.0 },
     energy_calc_update { 0.0 },
     traj_output_update { 0.0 },
     mpi_send_positions_update { 0.0 },
     mpi_send_forces_update { 0.0 },
     simulation_total { 0.0 },
     total_start { std::chrono::system_clock::now() } {}

    // Stop the total simulation timer and calculate the final timings.
    void finalize(void);

    // Start timers
    void start_cell_list_update(void);
    void start_force_update(void);
    void start_force_wall_update(void);
    void start_position_update(void);
    void start_velocity_update(void);
    void start_energy_calc_update(void);
    void start_traj_output_update(void);
    void start_mpi_send_positions_update(void);
    void start_mpi_send_forces_update(void);
    void start_mpi_clean_update(void);

    // Stop timers and add the duration from the start to the total.
    void stop_cell_list_update(void);
    void stop_force_update(void);
    void stop_force_wall_update(void);
    void stop_position_update(void);
    void stop_velocity_update(void);
    void stop_energy_calc_update(void);
    void stop_traj_output_update(void);
    void stop_mpi_send_positions_update(void);
    void stop_mpi_send_forces_update(void);
    void stop_mpi_clean_update(void);

    // Total time spent in different tasks.
    std::chrono::duration<double> cell_list_update,
                                  force_update,
                                  force_wall_update,
                                  position_update,
                                  velocity_update,
                                  energy_calc_update,
                                  traj_output_update,
                                  mpi_send_positions_update,
                                  mpi_send_forces_update,
                                  mpi_clean_update,
                                  simulation_total,
                                  rest;

private:
    ChronoTime cell_list_start,
               force_start,
               force_wall_start,
               position_start,
               velocity_start,
               energy_calc_start,
               traj_output_start,
               mpi_send_positions_start,
               mpi_send_forces_start,
               mpi_clean_start,
               total_start;
};

// Keep track of the non-dimensional energetics of the system during
// a simulation.
struct Energetics {
    std::vector<double> potential,
                        temperature,
                        times_per_calc;
    ChronoTime last_time;
};

// Describe the system.
void describe_system_config(const System& system, const ForceField& ff);

// Show the cell list distribution.
void show_atom_cell_list_distribution(const System& system);

// Print the benchmark data to stdout.
void print_benchmark(const Benchmark& bench, const uint64_t num_steps);

// Initiate calculation of system energetics.
void init_system_energetics(Energetics &energetics);

// Calculate the energetics of the system and append to the object vectors.
void calculate_system_energetics(Energetics& energetics,
                                 const System& system,
                                 const ForceField& ff);

void print_energetics(const Energetics& energy, const ForceField& ff);

// Calculate the mean velocity of the system.
RVec calc_mean_velocity(const System& system);

// Calculate the non-dimensional system temperature.
double calc_system_temperature(const System& system);

#endif // ANALYTICS_H
