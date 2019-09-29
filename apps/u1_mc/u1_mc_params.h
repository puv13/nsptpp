///! \file
/// Getopt paramters for u1_mc.cpp - Header

#pragma once

#include <getopt.h>
#include <string>

#ifdef _OPENMP
extern int n_threads; // Number of OMP threads to use
#endif

extern int hot_start;   // Hot or cold start for MC?
extern int tune_accept; // Tune MC acceptance rate?
extern int quiet;       // Generate less terminal output
extern int no_meas;     // Switch off measurements and only generate
                        // configs
extern double beta0;
extern std::string outfile;      // Name of outputfile
extern std::size_t measurements; // Number of measurenments
extern std::size_t swps;         // Number of MC sweeps between  measurements
extern std::size_t hit_param;    // Number of MC hits per site
extern std::size_t or_param;     // Number of OR sweeps per MC sweep
extern std::size_t thermal;      // Number of thermalisation sweeps
extern std::size_t Lt;           // Temporal extend of lattice
extern std::size_t Ls;           // Spacial extend of lattice

int parse_command_line_options(int argc, char *argv[]);
void print_params();
