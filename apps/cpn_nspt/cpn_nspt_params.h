///! \file
/// Getopt paramters for cpn_nspt.cpp - Header

#pragma once

#include <getopt.h>
#include <string>

#ifdef _OPENMP
extern int n_threads; // Number of OMP threads to use
#endif

extern int seq_init;       // Sequential (order by order) thermalisation?
extern int gaugefix;       // Switch to toggle stochastic gauge fixing
extern int zeromodes;      // Switch to toggle zero mode subtraction
extern int integrator;     // Type of integrator to use in Langevin?
extern int quiet;          // Generate less terminal output
extern int twisted;        // Determine boundary conditions
extern int debug;          // Switch on debug mode
extern int seed_only;      // Switch on seed-only mode
extern int use_seed_state; // Switch to decide if a stored RNG state
                           // is used when a seed lattice is loaded
extern int no_meas;        // Switch off measurements and only generate
                           // configs

extern double eps_t;                // Stochastic time step size
extern double eps_constr;           // Precision goal for the fulfilment of
                                    // constraints
extern std::string outfile;         // Name of output file
extern std::string seed_file;       // Name of seed file
extern std::size_t measurements;    // Number of measurenments
extern std::size_t constr_int;      // Number of steps between constraint
                                    // enforcements
extern std::size_t swps;            // Number of time steps between measurements
extern std::size_t thermal;         // Number of thermalisation time steps
extern std::size_t Lt;              // Temporal extend of lattice
extern std::size_t Ls;              // Spacial extend of lattice
extern std::string file_identifier; // A string added to the standard file name
                                    // to make it unique

int parse_command_line_options(int argc, char *argv[]);
void print_params();
