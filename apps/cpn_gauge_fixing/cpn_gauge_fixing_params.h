///! \file
/// Getopt paramters for cpn_gauge_fixing.cpp - Header

#pragma once

#include <getopt.h>
#include <string>

#ifdef _OPENMP
extern int n_threads; // Number of OMP threads to use
#endif
extern int quiet;            // Generate less terminal output
extern std::size_t num_gc;   // Number of gauge copies to consider for
                             // finding the global maximum
extern std::size_t max_iter; // Maximal number of iterations in the gauge
                             // fixing algorithm
extern double beta0;         // Choose a specific beta
extern double orparam;       // over-relaxation parameter
extern std::string infile;   // Name of input file
extern std::string outfile;  // Name of output file


int parse_command_line_options(int argc, char *argv[]);
void print_params();
