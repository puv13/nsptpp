///! \file
/// Getopt paramters for cpn_measurements.cpp - Header

#pragma once

#include <getopt.h>
#include <string>

#ifdef _OPENMP
extern int n_threads; // Number of OMP threads to use
#endif
extern int quiet;          // Generate less terminal output
extern int twisted;        // Switch on twisted BC
extern std::string infile; // Name of input file

int parse_command_line_options(int argc, char *argv[]);
void print_params();
