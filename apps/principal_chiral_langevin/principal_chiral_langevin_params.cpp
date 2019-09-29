///! \file
/// Getopt paramters for principal_chiral_langevin - Option parser

#include "principal_chiral_langevin_params.h"
#include <iomanip>
#include <iostream>

#ifdef _OPENMP
int n_threads = 0;
bool user_def = 0;
#endif

int hot_start = 1;
int quiet = 1;
int no_meas = 0;
int twisted = 0;
int debug = 0;
int integrator = 1;
int seed_only = 0;
int use_seed_state = 0;

std::size_t measurements = 100;
std::size_t swps = 5000;
std::size_t thermal = 200000;
std::size_t Lt = 10;
std::size_t Ls = 10;
double beta0 = 0.0;
double eps_t = 0.1;

std::string outfile = "";
std::string file_identifier = "";
std::string seed_file = "";
std::string intnames[] = { "Euler", "RK" };

///! Struct to define long options for getopt
static struct option long_options[] = {
  { "cold-start", no_argument, &hot_start, 0 },
  { "verbose", no_argument, &quiet, 0 },
  { "twisted", no_argument, &twisted, 1 },
  { "debug", no_argument, &debug, 1 },
  { "configs-only", no_argument, &no_meas, 1 },
  { "seed-only", no_argument, &seed_only, 1 },
  { "use-seed-state", no_argument, &use_seed_state, 1 },
  { "beta0", required_argument, NULL, 'b' },
  { "eps_t", required_argument, NULL, 'e' },
  { "measurements", required_argument, NULL, 'm' },
  { "sweeps", required_argument, NULL, 's' },
  { "thermalisation", required_argument, NULL, 't' },
  { "outfile", required_argument, NULL, 'f' },
  { "seed", required_argument, NULL, 'F' },
  { "lt", required_argument, NULL, 'T' },
  { "ls", required_argument, NULL, 'S' },
  { "num-threads", required_argument, NULL, 'n' },
  { "integrator", required_argument, NULL, 'i' },
  { "file-identifier", required_argument, NULL, 'x' },
  { 0, 0, 0, 0 }
};


static char short_option_list[] = "b:e:m:s:t:f:S:T:n:x:F:i:vw";

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
      case 'b':
        beta0 = static_cast<double>(std::stod(optarg));
        break;
      case 'e':
        eps_t = static_cast<double>(std::stod(optarg));
        break;
      case 'm':
        measurements = static_cast<std::size_t>(std::stoi(optarg));
        break;
      case 't':
        thermal = static_cast<std::size_t>(std::stoi(optarg));
        break;
      case 's':
        swps = static_cast<std::size_t>(std::stoi(optarg));
        break;
      case 'T':
        Lt = static_cast<std::size_t>(std::stoi(optarg));
        break;
      case 'S':
        Ls = static_cast<std::size_t>(std::stoi(optarg));
        break;
      case 'f':
        outfile = static_cast<std::string>(optarg);
        break;
      case 'F':
        seed_file = static_cast<std::string>(optarg);
        break;
      case 'i':
        integrator = static_cast<int>(std::stoi(optarg));
        break;
      case 'v':
        quiet = 0;
        break;
      case 'w':
        twisted = 1;
        break;
      case 'x':
        file_identifier = static_cast<std::string>(optarg);
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

  std::cout << "\tSimulation Parameters: " << std::endl;
  std::cout << "\t-------------------------------------------------------------"
               "---------"
            << std::endl;
  std::cout << std::left;
  std::cout << std::setw(f1_width)
            << "\tVerbose ouput:" << (quiet != 1 ? "ON" : "OFF") << std::endl;
  std::cout << std::setw(f1_width)
            << "\tHot Start: " << (hot_start == 1 ? "ON" : "OFF") << std::endl;
  std::cout << std::setw(f1_width) << "\tIntegrator: " << intnames[integrator]
            << std::endl;
  std::cout << std::setw(f1_width) << "\tTwisted BC: "
            << (twisted == 0 ? "OFF" : "FIRST SPATIAL COORD") << std::endl;
  std::cout << std::setw(f1_width)
            << "\tConfigs only:" << (no_meas != 1 ? "NO" : "YES") << std::endl;
  std::cout << std::setw(f1_width) << "\tThermalisation Sweeps: " << thermal
            << std::endl;
  std::cout << std::setw(f1_width) << "\tSweeps: " << swps << std::endl;
  std::cout << std::setw(f1_width)
            << (no_meas != 1 ? "\tMeasurements: " : "\tConfigurations: ")
            << measurements << std::endl;
  std::cout << std::setw(f1_width) << "\tStochastic time step: " << eps_t
            << std::endl;
  if (not outfile.empty()) {
    std::cout << std::setw(f1_width) << "\tOutput File: " << outfile
              << std::endl;
  }
  if (not seed_file.empty()) {
    std::cout << std::setw(f1_width) << "\tSeed File: " << seed_file
              << std::endl;
  }
  if (not file_identifier.empty()) {
    std::cout << std::setw(f1_width) << "\tFile ID: " << file_identifier
              << std::endl;
  }
  if (debug) {
    std::cout << std::setw(f1_width) << "\tDEBUG MODE: "
              << "ON" << std::endl;
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
