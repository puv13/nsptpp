///! \file
/// Read U1 configs from a file and fix them to a given gauge (at the moment
/// Landau gauge is the only gauge implemented).
///

#include "./gaugegroups/u1.h"
#include "../src/actions/compact_qed.h"
#include "../src/stat.h"
#include "./u1_gauge_fixing_params.h"
#include "io.h"
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <regex>
#include <string>
#include <vector>

// Global constants known at compile time
#ifdef SET_DIM_
const std::size_t dim = SET_DIM_;
#else
const std::size_t dim = 4;
#endif


// Allowed (relative) change in observables due to gauge fixing
double error_eps = 5 * 1.e-14;

int main(int argc, char *argv[])
{


  // ----------------------------------------------------------------------
  // Get parameters from command line
  // ----------------------------------------------------------------------
  parse_command_line_options(argc, argv);

  // ----------------------------------------------------------------------
  // Perform some sanity checks and print parameters
  // ----------------------------------------------------------------------
  if (not fileExists(infile)) {
    throw std::runtime_error("No valid input file specified.");
    return 0;
  }
  if (infile == outfile) {

    std::cout << "\n\tWARNING: Input and output file name should be different."
              << std::endl;
    outfile = "gf_" + outfile;
    std::cout << "\tWARNING: Renaming outfile to " << outfile << std::endl
              << std::endl;
  }
  if (outfile.empty()) {
    std::cout << "\n\tWARNING: No output file name given." << std::endl;
    outfile = "gf_" + infile;
    std::cout << "\tWARNING: Using default outfile name '" << outfile << "'\n"
              << std::endl;
  }
  print_params();

// Init OMP
#ifdef _OPENMP
  omp_set_dynamic(0);
  if (n_threads < 1) {
    n_threads = omp_get_max_threads();
  };
  omp_set_num_threads(n_threads);
#endif


  // ----------------------------------------------------------------------
  // Open input file (read-only)
  // ----------------------------------------------------------------------
  Hdf5File h5file(infile, true);

  // Get names of first level groups
  std::vector<std::string> obj_names = h5file.ObjectNames();

  // Init regex for group filtering and empty lattice
  std::regex beta_re("(beta)([[:digit:]]+\\.[[:digit:]]+)");
  LinkLattice<U1, dim> u1lat;

  // ----------------------------------------------------------------------
  // Loop over groups in infile and perform gauge fixing
  // ----------------------------------------------------------------------
  for (auto name : obj_names) {

    // std::cout << name << std::endl;

    // We want to match strings which fit the beta regex
    std::smatch matches;
    if (std::regex_search(name, matches, beta_re)) {

      // Ignore all other beta values if a specific beta0 was specified
      double beta = std::stod(matches[2].str());
      if (0. != beta0 && beta != beta0) {
        continue;
      }

      // Open group
      h5file.push(name);

      // Get all configs belonging to the current group
      auto cfg_names = h5file.ObjectNames();
      std::size_t num_cfg = cfg_names.size();
      std::size_t cnt_cfg = 0;

      // Status info
      if (not quiet) {
        std::cout.precision(4);
        std::cout << "\tWorking on beta " << std::fixed << std::setw(7) << beta
                  << ":" << std::endl;
      }

      double average_iter = 0;
      std::vector<double> v_quality;
      std::vector<double> v_orig_quality;
      for (auto cname : cfg_names) {

        ++cnt_cfg;
        // if ("cfg_0001" != cname){ continue; }

        // Read lattice
        h5file.read_lattice(u1lat, std::string(cname));

        // Read known value of some observables. This will be used to
        // cross-check the gauge fixing later on.
        double energy;
        h5file.readAttribute(energy, "energy_density", cname);

        // Init Action, observables and gauge fixing
        QEDAction<dim> u1_act(u1lat, beta);
        QEDGaugeFix<dim> u1_gf(u1lat, beta);

        // Compute 'gauge quality' of initial config and potential maximum
        auto LGQ = u1_gf.LandauGaugeQuality();
        v_orig_quality.push_back(LGQ);
        auto LG_MAX = static_cast<double>(u1lat.volume() * dim);

        // Fix to Landau Gauge
        std::size_t num_it = 0;
        if (LGQ > 1.e-9) {
          num_it = u1_gf.LandauGaugeDriver(num_gc, max_iter, orparam);
          average_iter += static_cast<double>(num_it);
        }
        v_quality.push_back(u1_gf.LandauGaugeQuality());

        // Compute some observables in the new gauge
        auto energy_obs = u1_act.energyDensity();

        // Cross check. Gauge transformations should not change the
        // value of gauge invariant observables.

        auto e_diff_rel = std::abs((energy - energy_obs) / energy);

        auto functional = u1_gf.LandauGaugeFunctional();
        auto quality = u1_gf.LandauGaugeQuality();

        if (not quiet) {
          std::cout.precision(4);
          std::cout << std::fixed << "\t\tconfiguration " << std::setw(4)
                    << std::right << cnt_cfg << "/" << std::left << std::setw(4)
                    << num_cfg << "\titer: " << std::right << std::setw(5)
                    << num_it << "\tBest Func.: " << std::setw(11)
                    << functional / LG_MAX << std::scientific
                    << "\tQuality: " << std::setw(11) << quality
                    << "\tOriginal: " << std::setw(11) << LGQ << std::endl;
        }

        // std::cout << "\tEnergy    Difference: " << e_diff_rel << std::endl;
        // std::cout << "Plaquette Difference: " << p_diff_rel << std::endl;

        if (e_diff_rel > error_eps) {
          std::cerr << "Error: Energy density changed during gauge fixing. "
                    << "Relative difference: " << std::scientific << e_diff_rel
                    << std::endl;
        }

        // Write result to outfile
        Hdf5File h5ofile(outfile);
        h5ofile.push(name);
        h5ofile.write_lattice(u1lat, cname);

        // Write some gauge fixing info
        h5ofile.addAttribute("Number of gauge copies", num_gc, cname);
        h5ofile.addAttribute("Maximum iterations", max_iter, cname);
        h5ofile.addAttribute("Relative gauge functional", functional / LG_MAX,
                             cname);
        h5ofile.addAttribute("Gauge quality", quality, cname);
        h5ofile.addAttribute("OR parameter", orparam, cname);
        h5ofile.addAttribute("energy_density", energy_obs, cname);
      }
      average_iter /= static_cast<double>(num_cfg * num_gc);
      if (not quiet) {
        std::cout.precision(2);
        std::cout
          << "\t\tDone! Average gauge fixing iterations per configuration: "
          << std::fixed << average_iter << std::endl;
        std::cout.precision(4);
        std::cout << "\t\t      Average gauge quality: " << std::scientific
                  << vec_mean(v_quality) << std::endl;
        std::cout << "\t\t      Average original gauge quality: "
                  << std::scientific << vec_mean(v_orig_quality) << std::endl;
      }

      h5file.pop();
    }
  }

  return 0;
}
