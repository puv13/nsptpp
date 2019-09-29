///! \file
/// Read CP(N) configs from a file and perform measurements
///

#include "./latticefields/cp.h"
#include "../src/actions/canonical_cpn.h"
#include "../src/stat.h"
#include "./cpn_measurements_params.h"
#include "./gaugegroups/u1.h"
#include "io.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <vector>


// Global constants known at compile time
const std::size_t dim = 2;

int main(int argc, char *argv[])
{

  // ----------------------------------------------------------------------
  // Init Random number generator
  // ----------------------------------------------------------------------
  std::random_device rd;
  std::uniform_int_distribution<size_t> udint;
  // std::size_t udint_max=udint.max();
  std::vector<size_t> seeds;
  for (auto i = 0ul; i < std::mt19937_64::state_size; ++i) {
    seeds.push_back(udint(rd));
  }

  std::seed_seq sseq(seeds.begin(), seeds.end());
  std::mt19937_64 generator(sseq);


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

  std::string suffix = "";
  if (0 != twisted) {
    suffix = suffix + "_twist";
  }

  Hdf5File h5file(infile, true);

  std::vector<std::string> obj_names = h5file.ObjectNames();
  FullLattice<CP<std::complex<double>, N>, U1, dim> cplat;

  bool init_rf = false;

  // ----------------------------------------------------------------------
  // Boundary Cond.
  // ----------------------------------------------------------------------
  // Boundary conditions
  // "Gauge field"
  std::array<double, dim> u1_bc_ar;
  u1_bc_ar.fill(0.0);
  BoundaryCondition<double, dim> u1_bc(u1_bc_ar, "PERIODIC U1");
  std::string cpn_str = "PERIODIC CPN";


  CPtwistPhase<std::complex<double>, N> no_twist(std::complex<double>(0.0));
  std::array<CPtwistPhase<std::complex<double>, N>, dim> tp_ar;
  tp_ar.fill(no_twist);

  if (twisted == 1) {

    CPtwistPhase<std::complex<double>, N> twist(std::complex<double>(0.0));

    for (auto i = 0LU; i < N; ++i) {
      twist[i] = static_cast<double>(i);
    }

    // Multiply all phases by 2pi/N
    twist = (std::atan(1.) * 8. / static_cast<double>(N)) * twist;

    // Twist in 1-direction
    tp_ar[1] = twist;
  }

  BoundaryCondition<CPtwistPhase<std::complex<double>, N>, dim> cpn_bc(tp_ar,
                                                                       cpn_str);


  // ----------------------------------------------------------------------
  // Loop over beta values
  // ----------------------------------------------------------------------
  for (auto name : obj_names) {
    // std::cout << name << std::endl;

    std::smatch matches;

    // Vectors to store results of measurements
    std::vector<double> v_energy;
    std::vector<std::vector<double>> v_ww_cor;
    std::vector<double> v_plaquette;
    std::vector<double> v_xi_G;
    std::vector<double> v_chi_m;
    std::vector<double> v_Q;
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
        h5file.read_lattice(cplat, std::string(cname));

        // Filename and header for output files
        // ----------------------------------------------------------------------
        std::stringstream outs;

        outs.str("");
        outs << "N" << std::defaultfloat << std::setfill('0') << std::setw(2)
             << N << "." << std::setfill('0') << std::setw(2)
             << cplat.dimensionsArray()[0] << "x" << std::setfill('0')
             << std::setw(2) << cplat.dimensionsArray()[1] << ".beta"
             << std::setw(5) << std::fixed << std::setprecision(3) << beta
             << suffix;

        std::string fname = outs.str() + "_obs.dat";
        std::string pfname = outs.str() + "_prp.dat";
        std::string wfname = outs.str() + "_wwc.dat";
        std::ofstream ofile;
        outs.str("");

        // Header for files
        // ----------------------------------------------------------------------
        if (first) {
          std::time_t result = std::time(NULL);
          outs << "\n\n# " << std::asctime(std::localtime(&result))
               << "# N=" << N << "\n"
               << "# Lt=" << cplat.dimensionsArray()[0] << "\n"
               << "# Lx=" << cplat.dimensionsArray()[1] << "\n"
               << "# beta=" << beta << "\n";


          // Scalar Obs
          ofile.open(fname, std::ios::app);
          ofile << outs.str() +
                     "# cfg_name\tEnergyDensity\tMeanPlaq\txi_G\tchi_m\tQ\n";
          ofile.close();

          // Props
          ofile.open(pfname, std::ios::app);
          ofile << outs.str() +
                     "# cfg_name\ttime\t<A0(t)A0(0)>\t<A1(t)A1(0)>\n";
          ofile.close();
          outs.str("");

          // WW Correlators
          ofile.open(wfname, std::ios::app);
          ofile << outs.str() + "# cfg_name\ttime\tG_W(t)\n";
          ofile.close();
          outs.str("");
        }

        // Measurements
        // ----------------------------------------------------------------------
        CanonicalCPNAction<std::complex<double>, N, dim> cp_act(cplat, beta);
        CanonicalCPNObservables<std::complex<double>, N, dim> cp_obs(cplat,
                                                                     beta);
        CanonicalCPNGaugeFix<std::complex<double>, N, dim> cp_gf(cplat, beta);

        v_energy.push_back(cp_act.energyDensity(cpn_bc));
        v_plaquette.push_back(cp_obs.meanPlaquette(u1_bc));
        v_xi_G.push_back(cp_obs.xi_G());
        v_chi_m.push_back(cp_obs.chi_m());
        v_Q.push_back(cp_obs.Q_top(cpn_bc));
        v_ww_cor.push_back(cp_obs.WW_correlators());


        outs << std::setfill(' ') << std::setprecision(6) << std::setw(10)
             << std::fixed << cname << std::scientific << "\t" << std::setw(10)
             << v_energy.back() << "\t" << std::setw(13) << v_plaquette.back()
             << "\t" << std::setw(13) << v_xi_G.back() << "\t" << std::setw(13)
             << v_chi_m.back() << "\t" << std::setprecision(2) << std::setw(6)
             << std::fixed << v_Q.back() << std::endl;

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
        auto GaugeFieldS = cp_obs.GaugeFieldSpatialSum();

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

        // Write WW correlators to file
        auto wwc = v_ww_cor.back();
        for (auto t = 0LU; t < wwc.size(); ++t) {
          outs.str("");
          outs << std::setfill(' ') << std::setprecision(6) << std::setw(10)
               << std::fixed << cname << std::scientific << "\t" << std::setw(4)
               << t << "\t" << std::setw(13) << wwc[t] << std::endl;

          ofile.open(wfname, std::ios::app);
          ofile << outs.str();
          ofile.close();
        }


        // Set bit to show that everything was initialised in first
        // loop
        first = false;
      }


      // for (auto it = At[0].begin();  it != At[0].end(); ++it)
      // {
      // 	std::cout << vec_mean((*it)) << "\t" ;
      // }


      // for (auto i = 0L; i< ((At[0])[0]).size(); ++i){

      // 	 auto aux0 =  (At[0][0])[i];
      // 	 auto aux1 =  (At[1][0])[i];

      // 	 for (auto j = 0L; j< At[0].size(); ++j)
      // 	 {
      // 	     //std::cout << "i: " << i << "\tj:" << j << std::endl;

      // 	     (At[0][j])[i] *= aux0;
      // 	     (At[1][j])[i] *= aux1;
      // 	 }
      // }


      // std::cout << "****\n";
      // for (auto it = At[1].begin();  it != At[1].end(); ++it)
      // {
      // 	std::cout << vec_mean((*it)) << "\t" ;
      // }


      // Do statistics
      auto mean_energy = vec_mean(v_energy);
      auto std_err_energy = vec_std_err(v_energy);
      auto mean_plaq = vec_mean(v_plaquette);
      auto std_err_plaq = vec_std_err(v_plaquette);
      auto mean_xi = vec_mean(v_xi_G);
      auto std_err_xi = vec_std_err(v_xi_G);
      auto mean_chi = vec_mean(v_chi_m);
      auto std_err_chi = vec_std_err(v_chi_m);
      auto mean_ww_cor = vec_vec_mean((v_ww_cor));
      auto std_err_ww_cor = vec_vec_std_err((v_ww_cor));
      auto mean_Q = vec_mean(v_Q);
      auto std_err_Q = vec_std_err(v_Q);
      auto chi_t = vec_exp2(v_Q) / static_cast<double>(cplat.volume());
      auto std_err_chi_t =
        vec_exp2_std_err(v_Q) / static_cast<double>(cplat.volume());

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
      std::string wresstr;
      std::ofstream resfile;
      std::ofstream presfile;
      std::ofstream wresfile;

      outs.str("");
      outs << "N" << std::setfill('0') << std::setw(2) << N << "."
           << std::setfill('0') << std::setw(2) << cplat.dimensionsArray()[0]
           << "x" << std::setfill('0') << std::setw(2)
           << cplat.dimensionsArray()[1];

      resstr = outs.str() + suffix + "_obs.dat";
      presstr = outs.str() + suffix + "_prp.dat";
      wresstr = outs.str() + suffix + "_wwc.dat";

      if (not init_rf) {
        std::time_t result = std::time(NULL);

        resfile.open(resstr, std::ios::app);
        resfile << "\n\n# " << std::asctime(std::localtime(&result))
                << "# beta\tEnergyDensity\tError\tMeanPlaq\tError\txi_G"
                   "\tError\tchi_m\tError\tQ\tError\tAccCPN\tAccU1\n";
        resfile.close();

        presfile.open(presstr, std::ios::app);
        presfile
          << "\n\n# " << std::asctime(std::localtime(&result))
          << "# beta   \ttime\t<A0(t)A0(0)>\tError\t<A1(t)A1(0)>\tError\n";
        presfile.close();

        wresfile.open(wresstr, std::ios::app);
        wresfile << "\n\n# " << std::asctime(std::localtime(&result))
                 << "# beta   \ttime\tG_w(t)\tError\n";
        wresfile.close();

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
           << std::setw(13) << mean_xi << "\t" << std::setw(13) << std_err_xi
           << "\t" << std::setw(13) << mean_chi << "\t" << std::setw(13)
           << std_err_chi << "\t" << std::setw(13) << mean_Q << "\t"
           << std::setw(13) << std_err_Q << "\t" << std::setw(13) << chi_t
           << "\t" << std::setw(13) << std_err_chi_t << "\t" << std::endl;

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

      // WW Corr
      wresfile.open(wresstr, std::ios::app);
      outs.str("");
      for (auto t = 0LU; t < mean_ww_cor.size(); ++t) {
        outs << std::setfill(' ') << std::fixed << std::setprecision(4)
             << std::setw(8) << beta << "\t" << std::scientific
             << std::setprecision(6) << std::setw(04) << t << "\t"
             << std::setw(13) << mean_ww_cor[t] << "\t" << std::setw(13)
             << std_err_ww_cor[t] << std::endl;
      }

      outs << std::endl << std::endl;
      wresfile << outs.str();
      wresfile.close();
      outs.str("");


      // return 0.;
      h5file.pop();
    }
  }

  return 0;
}
