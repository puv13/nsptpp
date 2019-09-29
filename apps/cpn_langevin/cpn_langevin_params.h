///! \file
/// Getopt paramters for cpn_langevin.cpp - Header

#pragma once

#include <getopt.h>
#include <string>

#ifdef _OPENMP
extern int n_threads; // Number of OMP threads to use
#endif

extern int hot_start;  // Hot or cold start for MC?
extern int integrator; // Type of integrator to use in Langevin?
extern int quiet;      // Generate less terminal output
extern int twisted;    // Determine boundary conditions
extern int debug;      // Switch on debug mode
extern int no_meas;    // Switch off measurements and only generate
                       // configs
extern double beta0;
extern double eps_t;             // Stochastic time step size
extern std::string outfile;      // Name of outputfile
extern std::size_t measurements; // Number of measurenments
extern std::size_t swps;         // Number of time steps between measurements
extern std::size_t thermal;      // Number of thermalisation time steps
extern std::size_t Lt;           // Temporal extend of lattice
extern std::size_t Ls;           // Spacial extend of lattice

int parse_command_line_options(int argc, char *argv[]);
void print_params();
