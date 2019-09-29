///! \file
/// Read CP(N) configs from a file and fix them to a given gauge (at the moment
/// Landau gauge is the only gauge implemented).
///

#include "./latticefields/cp.h"
#include "../src/actions/canonical_cpn.h"
#include "./cpn_gauge_fixing_params.h"
#include "./gaugegroups/u1.h"
#include "io.h"
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <regex>
#include <string>
#include <vector>

// Global constants known at compile time
const std::size_t dim = 2;

// Allowed (relative) change in observables due to gauge fixing
double error_eps = 5 * N * 1.e-14;

int main(int argc, char *argv[])
{


  // ----------------------------------------------------------------------
  // Get parameters from command line
  // ----------------------------------------------------------------------
  parse_command_line_options(argc, argv);

  // ----------------------------------------------------------------------
  // Perform some sanity checks and print parameters
  // ----------------------------------------------------------------------
  std::cout << "\tThis program was compiled for CP(N-1) with N=" << N
            << std::endl
            << std::endl;
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
  FullLattice<CP<std::complex<double>, N>, U1, dim> cplat;

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
      for (auto cname : cfg_names) {

        ++cnt_cfg;
        // if ("cfg_0001" != cname){ continue; }

        // Read lattice
        h5file.read_lattice(cplat, std::string(cname));

        // Read known value of some observables. This will be used to
        // cross-check the gauge fixing later on.
        double energy;
        h5file.readAttribute(energy, "energy_density", cname);
        double plaq;
        h5file.readAttribute(plaq, "mean_plaquette", cname);

        // Init Action, observables and gauge fixing
        CanonicalCPNAction<std::complex<double>, N, dim> cp_act(cplat, beta);
        CanonicalCPNObservables<std::complex<double>, N, dim> cp_obs(cplat,
                                                                     beta);
        CanonicalCPNGaugeFix<std::complex<double>, N, dim> cp_gf(cplat, beta);

        // Compute 'gauge quality' of initial config and potential maximum
        auto LGQ = cp_gf.LandauGaugeQuality();
        auto LG_MAX = static_cast<double>(cplat.volume() * dim);

        // Fix to Landau Gauge
        std::size_t num_it = cp_gf.LandauGaugeDriver(num_gc, max_iter, orparam);
        average_iter += static_cast<double>(num_it);

        // Compute some observables in the new gauge
        auto energy_obs = cp_act.energyDensity();
        auto plaq_obs = cp_obs.meanPlaquette();

        // Cross check. Gauge transformations should not change the
        // value of gauge invariant observables.

        auto e_diff_rel = std::abs((energy - energy_obs) / energy);
        auto p_diff_rel = std::abs((plaq - plaq_obs) / plaq);

        auto functional = cp_gf.LandauGaugeFunctional();
        auto quality = cp_gf.LandauGaugeQuality();

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

        // std::cout << "Energy    Difference: " << e_diff_rel << std::endl;
        // std::cout << "Plaquette Difference: " << p_diff_rel << std::endl;

        if (e_diff_rel > error_eps) {
          std::cerr << "Error: Energy density changed during gauge fixing. "
                    << "Relative difference: " << std::scientific << e_diff_rel
                    << std::endl;
        }
        if (p_diff_rel > error_eps) {
          std::cerr << "Error: Mean plaquette changed during gauge fixing. "
                    << "Relative difference: " << std::scientific << p_diff_rel
                    << std::endl;
        }

        // Write result to outfile
        Hdf5File h5ofile(outfile);
        h5ofile.push(name);
        h5ofile.write_lattice(cplat, cname);

        // Write some gauge fixing info
        h5ofile.addAttribute("Number of gauge copies", num_gc, cname);
        h5ofile.addAttribute("Maximum iterations", max_iter, cname);
        h5ofile.addAttribute("Relative gauge functional", functional / LG_MAX,
                             cname);
        h5ofile.addAttribute("Gauge quality", quality, cname);
        h5ofile.addAttribute("OR parameter", orparam, cname);
      }
      average_iter /= static_cast<double>(num_cfg * num_gc);
      if (not quiet) {
        std::cout.precision(2);
        std::cout
          << "\t\tDone! Average gauge fixing iterations per configuration: "
          << std::fixed << average_iter << std::endl;
      }

      h5file.pop();
    }
  }

  return 0;
}

// for (size_t i = 0; i < matches.size(); ++i) {
// 	 std::cout << i << ": '" << matches[i].str() << "'\n";
// }
