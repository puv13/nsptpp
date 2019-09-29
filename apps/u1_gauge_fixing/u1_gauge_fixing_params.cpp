///! \file
/// Getopt paramters for u1_gauge_fixing.cpp - Option parser
#include "./u1_gauge_fixing_params.h"
#include <iomanip>
#include <iostream>


#ifdef _OPENMP
int n_threads = 0;
bool user_def = 0;
#endif
int quiet = 1;
std::size_t num_gc = 1;
std::size_t max_iter = 1000;
double beta0 = 0.0;
double orparam = 1.0;
std::string infile = "";
std::string outfile = "";


///! Struct to define long options for getopt
static struct option long_options[] = {
  { "verbose", no_argument, &quiet, 0 },
  { "infile", required_argument, NULL, 'f' },
  { "outfile", required_argument, NULL, 'o' },
  { "beta", required_argument, NULL, 'b' },
  { "or-param", required_argument, NULL, 'r' },
  { "num-gc", required_argument, NULL, 'c' },
  { "max_iter", required_argument, NULL, 'i' },
  { "num-threads", required_argument, NULL, 'n' },
  { 0, 0, 0, 0 }
};


static char short_option_list[] = "m:s:t:f:o:b:r:c:i:n:v";

int parse_command_line_options(int argc, char **argv)
{
  int c, option_index;

  while (1) {
    c = getopt_long(argc, argv, short_option_list, long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
      case 0:
        break;
      case 'f':
        infile = static_cast<std::string>(optarg);
        break;
      case 'o':
        outfile = static_cast<std::string>(optarg);
        break;
      case 'b':
        beta0 = static_cast<double>(std::stod(optarg));
        break;
      case 'r':
        orparam = static_cast<double>(std::stod(optarg));
        break;
      case 'c':
        num_gc = static_cast<std::size_t>(std::stoi(optarg));
        break;
      case 'i':
        max_iter = static_cast<std::size_t>(std::stoi(optarg));
        break;
#ifdef _OPENMP
      case 'n':
        n_threads = static_cast<int>(std::stoi(optarg));
        if (n_threads >= 0) {
          user_def = true;
        }
        else {
          n_threads *= -1;
        }
        break;
#endif
      case 'v':
        quiet = 0;
        break;

      default:
        break;
    }
  }

  return 0;
}


void print_params()
{

  const int f1_width = 28;

  std::cout << "\tParameters: " << std::endl;
  std::cout << "\t-------------------------------------------------------------"
               "---------"
            << std::endl;
  std::cout << std::left;
  std::cout << std::setw(f1_width)
            << "\tVerbose ouput:" << (quiet != 1 ? "ON" : "OFF") << std::endl;
  if (not infile.empty()) {
    std::cout << std::setw(f1_width) << "\tInput File: " << infile << std::endl;
  }
  if (not outfile.empty()) {
    std::cout << std::setw(f1_width) << "\tOutput File: " << outfile
              << std::endl;
  }
  if (0.0 != beta0) {
    std::cout << std::setw(f1_width) << "\tBeta: " << beta0 << std::endl;
  }
  std::cout << std::setw(f1_width) << "\tGauge Copies: " << num_gc << std::endl;
  std::cout << std::setw(f1_width) << "\tMax. Iterations: " << max_iter
            << std::endl;
  std::cout << std::setw(f1_width) << "\tOR parameter: " << orparam
            << std::endl;
#ifdef _OPENMP
  if (user_def) {
    std::cout << std::setw(f1_width) << "\tOMP threads:" << n_threads
              << std::endl;
  }
  else {
    std::cout << std::setw(f1_width) << "\tOMP threads:"
              << "AUTOMATIC" << std::endl;
  }
#endif
  std::cout << "\t-------------------------------------------------------------"
               "---------\n"
            << std::endl;
}
