///! \file
/// Getopt paramters for u1_measurements.cpp - Header

#pragma once

#include <getopt.h>
#include <string>

#ifdef _OPENMP
extern int n_threads; // Number of OMP threads to use
#endif
extern int quiet;    // Generate less terminal output
extern double beta0; // Choose a specific beta

extern std::string infile; // Name of input file


int parse_command_line_options(int argc, char *argv[]);
void print_params();
