///! \file
/// Read U1 configs from a file and perform measurements
///

#include "./gaugegroups/u1.h"
#include "../src/actions/compact_qed.h"
#include "../src/stat.h"
#include "./u1_measurements_params.h"
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

  std::regex beta_re("(beta)([[:digit:]]+\\.[[:digit:]]+)");

  // ----------------------------------------------------------------------
  // Get parameters from command line
  // ----------------------------------------------------------------------
  parse_command_line_options(argc, argv);
  print_params();

// Init OMP
#ifdef _OPENMP
  omp_set_dynamic(0);
  if (n_threads < 1) {
    n_threads = omp_get_max_threads();
  };
  omp_set_num_threads(n_threads);
#endif

  Hdf5File h5file(infile, true);

  std::vector<std::string> obj_names = h5file.ObjectNames();
  LinkLattice<U1, dim> u1_lat;

  bool init_rf = false;

  // ----------------------------------------------------------------------
  // Boundary Cond.
  // ----------------------------------------------------------------------
  // Boundary conditions
  // "Gauge field"
  std::array<double, dim> u1_bc_ar;
  u1_bc_ar.fill(0.0);
  BoundaryCondition<double, dim> u1_bc(u1_bc_ar, "PERIODIC U1");


  // ----------------------------------------------------------------------
  // Loop over beta values
  // ----------------------------------------------------------------------
  for (auto name : obj_names) {
    // std::cout << name << std::endl;

    std::smatch matches;

    // Vectors to store results of measurements
    std::vector<double> v_energy;
    std::vector<double> v_plaquette;
    // Storage for gauge props
    std::array<std::vector<std::vector<double>>, dim> At;

    if (std::regex_search(name, matches, beta_re)) {
      // Open group
      h5file.push(name);
      double beta = std::stof(matches[2].str());

      // Get content
      auto cfg_names = h5file.ObjectNames();


      bool first = true;
      for (auto cname : cfg_names) {

        // if ("cfg_0001" != cname){ continue; }

        // Read Lattice
        h5file.read_lattice(u1_lat, std::string(cname));


        // Filename and header for output files
        // ----------------------------------------------------------------------
        std::stringstream outs;

        outs.str("");
        outs << "cQED" << dim << "." << std::setfill('0') << std::setw(2)
             << u1_lat.dimensionsArray()[0] << "x" << std::setfill('0')
             << std::setw(2) << u1_lat.dimensionsArray()[1] << ".beta"
             << std::setw(5) << std::fixed << std::setprecision(3) << beta;

        std::string fname = outs.str() + "_obs.dat";
        std::string pfname = outs.str() + "_prp.dat";

        std::ofstream ofile;
        outs.str("");

        // Header for files
        // ----------------------------------------------------------------------
        if (first) {
          std::time_t result = std::time(NULL);
          outs << "\n\n# " << std::asctime(std::localtime(&result))
               << "# Lt=" << u1_lat.dimensionsArray()[0] << "\n"
               << "# Lx=" << u1_lat.dimensionsArray()[1] << "\n"
               << "# beta=" << beta << "\n";


          // Scalar Obs
          ofile.open(fname, std::ios::app);
          ofile << outs.str() + "# cfg_name\tEnergyDensity\tMeanPlaq\n";
          ofile.close();

          // Props
          ofile.open(pfname, std::ios::app);
          ofile << outs.str() +
                     "# cfg_name\ttime\t<A0(t)A0(0)>\t<A1(t)A1(0)>\n";
          ofile.close();
          outs.str("");
        }

        // Measurements
        // ----------------------------------------------------------------------
        QEDAction<dim> u1_act(u1_lat, beta);

        v_energy.push_back(u1_act.energyDensity(u1_bc));
        v_plaquette.push_back(u1_act.meanPlaquette(u1_bc));


        outs << std::setfill(' ') << std::setprecision(6) << std::setw(10)
             << std::fixed << cname << std::scientific << "\t" << std::setw(10)
             << v_energy.back() << "\t" << std::setw(13) << v_plaquette.back()
             << std::endl;

        ofile.open(fname, std::ios::app);
        ofile << outs.str();
        ofile.close();

        if (!quiet) {
          std::cout << outs.str();
        }
        outs.str("");

        // Compute propagators
        // ----------------------------------------------------------------------

        // Get A0(t) and A1(t)
        auto GaugeFieldS = u1_act.GaugeFieldSpatialSum();

        std::size_t cnt = 0;

        // Propagators
        for (auto it = GaugeFieldS.begin(); it != GaugeFieldS.end(); ++it) {
          ++cnt;

          // Store A0(0)  and A1(0)
          double aux[dim];
          if (GaugeFieldS.begin() == it) {
            for (auto i = 0LU; i < dim; ++i) {
              aux[i] = (*it)[i];
            }
          }

          // Compute propagators from summed gauge fields and store
          // results
          for (auto i = 0UL; i < dim; ++i) {

            if (first) {
              std::vector<double> vec = { (*it)[i] * aux[i] };
              At[i].push_back(vec);
              // std::cout << "First: " << At[i].size() << "\t" ;
            }
            else {
              // std::cout << At[i].size() << "/" << cnt << std::endl;
              ((At[i])[cnt - 1]).push_back((*it)[i] * aux[i]);
            }
          }
        }

        // Write props to file
        for (auto t = 0LU; t < At[0].size(); ++t) {
          outs.str("");
          outs << std::setfill(' ') << std::setprecision(6) << std::setw(10)
               << std::fixed << cname << std::scientific << "\t" << std::setw(4)
               << t << "\t" << std::setw(13) << (At[0])[t].back() << "\t"
               << std::setw(13) << (At[1])[t].back() << std::endl;

          ofile.open(pfname, std::ios::app);
          ofile << outs.str();
          ofile.close();
        }

        // Set bit to show that everything was initialised in first
        // loop
        first = false;
      }


      // Do statistics
      auto mean_energy = vec_mean(v_energy);
      auto std_err_energy = vec_std_err(v_energy);
      auto mean_plaq = vec_mean(v_plaquette);
      auto std_err_plaq = vec_std_err(v_plaquette);

      std::array<std::vector<double>, dim> mean_At;
      std::array<std::vector<double>, dim> std_err_At;

      // Statistics for propagators
      for (auto t = 0LU; t < At[0].size(); ++t) {
        for (auto i = 0LU; i < dim; ++i) {
          mean_At[i].push_back(vec_mean((At[i])[t]));
          // std::cout << "i/t: " << i << "/" << t << "\t" << vec_mean(
          // (At[i])[t]) << "\n";
          std_err_At[i].push_back(vec_std_err((At[i])[t]));
        }
        // std::cout << vec_mean((*it)) << "\t" ;
      }


      // Write results to file
      //----------------------------------------------------------------------
      std::stringstream outs;
      std::string resstr;
      std::string presstr;
      std::ofstream resfile;
      std::ofstream presfile;

      outs.str("");
      outs << "cQED" << dim << "." << std::setfill('0') << std::setw(2)
           << u1_lat.dimensionsArray()[0] << "x" << std::setfill('0')
           << std::setw(2) << u1_lat.dimensionsArray()[1];

      resstr = outs.str() + "_obs.dat";
      presstr = outs.str() + "_prp.dat";

      if (not init_rf) {
        std::time_t result = std::time(NULL);

        resfile.open(resstr, std::ios::app);
        resfile << "\n\n# " << std::asctime(std::localtime(&result))
                << "# beta\tEnergyDensity\tError\tMeanPlaq\tError\n";
        resfile.close();

        presfile.open(presstr, std::ios::app);
        presfile
          << "\n\n# " << std::asctime(std::localtime(&result))
          << "# beta   \ttime\t<A0(t)A0(0)>\tError\t<Ai(t)Ai(0)>\tError\n";
        presfile.close();

        init_rf = true;
      }

      // Scalar observables
      resfile.open(resstr, std::ios::app);

      outs.str("");
      outs << std::setfill(' ') << std::fixed << std::setprecision(4)
           << std::setw(8) << beta << "\t" << std::scientific
           << std::setprecision(6) << std::setw(13) << mean_energy << "\t"
           << std::setw(13) << std_err_energy << "\t" << std::setw(13)
           << mean_plaq << "\t" << std::setw(13) << std_err_plaq << "\t"
           << std::endl;

      resfile << outs.str();
      resfile.close();
      outs.str("");

      // Propagators
      presfile.open(presstr, std::ios::app);
      outs.str("");
      for (auto t = 0LU; t < mean_At[0].size(); ++t) {
        outs << std::setfill(' ') << std::fixed << std::setprecision(4)
             << std::setw(8) << beta << "\t" << std::scientific
             << std::setprecision(6) << std::setw(04) << t << "\t";

        for (auto i = 0LU; i < dim; ++i) {
          outs << std::setw(13) << mean_At[i][t] << "\t" << std::setw(13)
               << std_err_At[i][t] << "\t";
        }
        outs << std::endl;
      }

      outs << std::endl << std::endl;
      presfile << outs.str();
      presfile.close();
      outs.str("");


      // return 0.;
      h5file.pop();
    }
  }
  return 0;
}
