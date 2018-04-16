#include <chrono>

#include "conf.h"
#include "params.h"

#ifndef ANALYTICS_H
#define ANALYTICS_H

class Benchmark {
public:
    // Initialize all timers to 0 and start the total simulation timer.
    Benchmark(void)
    :cell_list_update { 0.0 },
     force_update { 0.0 },
     position_update { 0.0 },
     velocity_update { 0.0 },
     simulation_total { 0.0 },
     total_start { std::chrono::system_clock::now() } {}

    // Stop the total simulation timer and calculate the final timings.
    void finalize(void);

    // Start timers
    void start_cell_list_update(void);
    void start_force_update(void);
    void start_position_update(void);
    void start_velocity_update(void);

    // Stop timers and add the duration from the start to the total.
    void stop_cell_list_update(void);
    void stop_force_update(void);
    void stop_position_update(void);
    void stop_velocity_update(void);

    // Total time spent in different tasks.
    std::chrono::duration<double> cell_list_update,
                                  force_update,
                                  position_update,
                                  velocity_update,
                                  simulation_total,
                                  rest;

private:
    decltype(std::chrono::system_clock::now()) cell_list_start,
                                               force_start,
                                               position_start,
                                               velocity_start,
                                               total_start;
};

// Keep track of the energetics of the system during a simulation.
struct Energetics {
    std::vector<double> potential,
                        temperature;
};

void describe_system_config(const System& system);
void show_atom_cell_list_distribution(const System& system);
void print_benchmark(const Benchmark& bench);
void calculate_system_energetics(Energetics& energetics,
                                 const System& system,
                                 const ForceField& ff);

#endif // ANALYTICS_H
