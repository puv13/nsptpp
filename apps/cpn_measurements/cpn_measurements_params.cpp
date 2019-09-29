///! \file
/// Getopt paramters for cpn_measurements.cpp - Option parser


#include "./cpn_measurements_params.h"
#include <iomanip>
#include <iostream>

#ifdef _OPENMP
int n_threads = 0;
bool user_def = 0;
#endif
int quiet = 1;
int twisted = 0;
std::string infile = "";


///! Struct to define long options for getopt
static struct option long_options[] = {
  { "verbose", no_argument, &quiet, 0 },
  { "twisted", no_argument, &twisted, 1 },
  { "infile", required_argument, NULL, 'f' },
  { "num-threads", required_argument, NULL, 'n' },
  { 0, 0, 0, 0 }
};


static char short_option_list[] = "m:s:t:f:n:vw";

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
      case 'w':
        twisted = 1;
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
  std::cout << std::setw(f1_width) << "\tTwisted BC: "
            << (twisted == 0 ? "OFF" : "FIRST SPATIAL COORD") << std::endl;

  if (not infile.empty()) {
    std::cout << std::setw(f1_width) << "\tInput File: " << infile << std::endl;
  }
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
