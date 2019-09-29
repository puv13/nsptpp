///! \file
/// An implementation of the Langevin algorithm for the CP(N-1) model
///

#ifdef BLAS_AVAIL
#include "cpn_langevin.h"
#include "io.h"
#include "randgen.h"
#include <omp.h>


int main(int argc, char *argv[])
{

  // ----------------------------------------------------------------------
  // Say Hi
  // ----------------------------------------------------------------------

  std::cout << "\nRunning Langevin for the CP(N-1) model ..." << std::endl;
  std::cout << "________________________________________________"
            << "______________________________\n"
            << std::endl;

  // ----------------------------------------------------------------------
  // Get Langevin  parameters from command line options
  // ----------------------------------------------------------------------

  parse_command_line_options(argc, argv);
  print_params();

// ----------------------------------------------------------------------
// Init OMP
// ----------------------------------------------------------------------
#ifdef _OPENMP
  omp_set_dynamic(0);
  if (n_threads < 1) {
    n_threads = omp_get_max_threads();
  };
  omp_set_num_threads(n_threads);
#endif

  // // ----------------------------------------------------------------------
  // // Sanity check
  // // ----------------------------------------------------------------------
  // if (outfile.empty() && no_meas)
  // {
  // 	 std::cout <<
  // "**********************************************************************\n";
  // 	 std::cout << "You don't want to do measurements and "
  // 		   << "you didn't specify a file\nto safe the "
  // 		   << "configurations.\n"
  // 		   << "Why even bother?\n"
  // 		   << "... I am out ... " << std::endl;
  // 	 std::cout <<
  // "**********************************************************************\n";

  // 	 return 1;
  // }

  // ----------------------------------------------------------------------
  // Define betas by hand
  // ----------------------------------------------------------------------

  std::vector<double> scan_range;

  if (beta0 <= 0.0) {

    scan_range = std::vector<double>({
      //.8, .85,.9,.95
      0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
      1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2,
      2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0
      // 	 5.0 , 10. ,20.,30.,40.
    });
  }
  else {
    scan_range.push_back(beta0);
  }


  // ----------------------------------------------------------------------
  // Print Lattice setup
  // ----------------------------------------------------------------------
  std::cout << "\tLattice Parameters: " << std::endl;
  std::cout << "\t-------------------------------------------------------------"
               "---------"
            << std::endl;
  std::cout << std::left;
  std::cout << std::setw(28) << "\tLt:" << Lt << std::endl;
  std::cout << std::setw(28) << "\tLs:" << Ls << std::endl;
  std::cout << std::setw(28) << "\tN:" << N << std::endl;
  std::cout << std::setw(28) << "\tbeta:";
  if (1UL == scan_range.size()) {
    std::cout << scan_range.front() << std::endl;
  }
  else {
    int cnt = 0;
    for (auto beta_it = scan_range.begin(); beta_it != scan_range.end();
         ++beta_it) {
      beta_it != scan_range.begin() ? (std::cout << "; ") : (std::cout << "{");
      cnt++;
      if (0 == cnt % 5) {
        std::cout << std::endl << std::setw(36) << "";
      }
      std::cout << std::fixed << std::setprecision(2) << std::setw(5)
                << std::right << *beta_it;
    }
    std::cout << "}" << std::endl;
  }
  std::cout << "\t-------------------------------------------------------------"
               "---------\n"
            << std::endl;


  // ----------------------------------------------------------------------
  // Init Random number generator
  // ----------------------------------------------------------------------
  std::random_device rd;
  std::uniform_int_distribution<size_t> udint;
  std::vector<size_t> seeds;
  if (debug) {
    for (auto i = 0ul; i < std::mt19937_64::state_size; ++i) {
      seeds.push_back(i);
    }
  }
  else {
    for (auto i = 0ul; i < std::mt19937_64::state_size; ++i) {
      seeds.push_back(udint(rd));
    }
  }

  std::seed_seq sseq(seeds.begin(), seeds.end());
  std::mt19937_64 generator(sseq);

  // ----------------------------------------------------------------------
  // Init Boundaries
  // ----------------------------------------------------------------------

  // "Gauge field"
  std::array<double, dim> u1_bc_ar;
  u1_bc_ar.fill(1.0);
  BoundaryCondition<double, dim> u1_bc(u1_bc_ar, "PERIODIC U1");
  std::string cpn_str = "PERIODIC CPN";


  // CPN
  CPtwistPhase<std::complex<double>, N> no_twist(std::complex<double>(1.0));
  std::array<CPtwistPhase<std::complex<double>, N>, dim> tp_ar;
  tp_ar.fill(no_twist);

  if (twisted == 1) {

    CPtwistPhase<std::complex<double>, N> twist(std::complex<double>(0.0));

    auto ii = std::complex<double>(0., 1.);
    for (auto i = 0LU; i < N; ++i) {
      twist[i] = std::exp(-ii * (std::atan(1.) * 8. / static_cast<double>(N)) *
                          static_cast<double>(i));
    }


    tp_ar[1] = twist;
  }

  BoundaryCondition<CPtwistPhase<std::complex<double>, N>, dim> cpn_bc(tp_ar,
                                                                       cpn_str);


  // ----------------------------------------------------------------------
  // Resfile header
  // ----------------------------------------------------------------------
  std::string start_cond = (hot_start) ? "" : "_cs";
  std::string integr_str = (integrator > 0) ? "" : "_euler";
  if (0 != twisted) {
    start_cond = start_cond + "_twist";
  }
  std::ofstream resfile;
  std::ostringstream epsstr;
  epsstr << std::setprecision(4) << std::scientific << eps_t;

  std::stringstream outs;
  outs << "N" << std::setfill('0') << std::setw(2) << N << "."
       << std::setfill('0') << std::setw(2) << Lt << "x" << std::setfill('0')
       << std::setw(2) << Ls << "_eps" << epsstr.str() << start_cond
       << integr_str << ".dat";
  auto resstr = outs.str();
  outs.str("");
  resfile.open(resstr, std::ios::app);
  resfile << "\n\n# "
             "beta\tEnergyDensity\tError\tMeanPlaq\tError\txi_G\tError\tchi_"
             "m\tError\tQ\tError\tchi_t\tError\n";
  resfile.close();


  // ----------------------------------------------------------------------
  // Loop over beta values
  // ----------------------------------------------------------------------
  std::cout << "Looping over beta ..." << std::endl;
  std::cout << "_______________________________________________________________"
               "_______________\n"
            << std::endl;


  for (auto beta : scan_range) {

    std::cout << "\n\tWorking on beta=" << beta << std::endl;
    std::cout << "\t..........................................................."
                 "..........."
              << std::endl;
    auto begin = clck::now();

    // ----------------------------------------------------------------------
    // Langevin Driver
    // ----------------------------------------------------------------------
    CPN cp_init(1.0, true);
    U1 u1_init(-1.0);

    ranlux24_normal<double> rlx;

    // This is DEBUG stuff!! Same noise every time ...
    if (debug) {
      rlx.set_generator_state("1857670 9122615 9091361 7165633 6622926 5093178 "
                              "10658136 14316459 2246149 476568 820330 "
                              "11401107 5667087 14569279 13287542 9551677 "
                              "2350107 409682 3168448 8334688 10338889 "
                              "15771630 12488960 15078927 0 0 0");
      rlx.set_distribution_state(
        "0.00000000000000000e+00 1.00000000000000000e+00 0");
    }


    std::array<size_t, dim> dimsar;
    dimsar.fill(Ls);
    dimsar[0] = Lt;
    CPlattice lattice(dimsar, cp_init, u1_init);
    CanonicalCPNAction<std::complex<double>, N, dim> cp_act(lattice, beta);
    CanonicalCPNObservables<std::complex<double>, N, dim> cp_obs(lattice, beta);

    NoExpansion::CPForce<N, dim, CPtwistPhase<std::complex<double>, N>, double>
      Force(cpn_bc, u1_bc);


    IntegratorsNoExp::EulerInt<
      CPlattice, typename std::vector<CPForceStruct<N, dim>>, decltype(Force),
      randomCPdf<ranlux24_normal<double>>, ranlux24_normal<double>>
      euler;
    euler.setBeta(beta);
    euler.setStepSize(eps_t);

    IntegratorsNoExp::RKInt<
      CPlattice, typename std::vector<CPForceStruct<N, dim>>,
      CPCasimirStruct<N>, decltype(Force), randomCPdf<ranlux24_normal<double>>,
      ranlux24_normal<double>>
      rkint;
    rkint.setBeta(beta);
    rkint.setStepSize(eps_t);

    // Hot or cold start?
    if (hot_start) {
      canonical_CP_random_init(lattice, generator);
    }


    // Thermalisation
    // ----------------------------------------------------------------------
    std::cout << "\tThermalisation ... " << std::endl;
    auto start = clck::now();

    // Show progress of thermalisation
    const size_t prog_steps = 100;
    progress_bar pb(static_cast<double>(prog_steps), 56, std::string("\t\t"));
    for (auto t = 0LU; t < prog_steps; ++t) {
      pb.progress(static_cast<double>(t));
      size_t th = static_cast<size_t>(std::ceil(
        static_cast<double>(thermal) / static_cast<double>(prog_steps)));
      for (auto i = 0LU; i < th; ++i) {
        if (integrator > 0) {
          rkint.update(Force, lattice, rlx);
        }
        else {
          euler.update(Force, lattice, rlx);
        }
      }
    }
    pb.progress(static_cast<double>(prog_steps));
    std::cout << std::endl;

    auto end = clck::now();
    auto time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "\t\t... done\n" << std::endl;


    // Measurements & Config Generation
    // ----------------------------------------------------------------------

    // Observables
    // Vectors to store results of measurements
    std::vector<double> v_energy;
    std::vector<double> v_plaquette;
    std::vector<double> v_xi_G;
    std::vector<double> v_chi_m;
    std::vector<double> v_Q;

    int m_length;
    std::ofstream ofile;
    std::string fname;

    if (not no_meas) {

      std::cout << "\tMeasurements ... " << std::endl;


      // Filename and header for file output
      outs.str("");
      outs << "N" << std::defaultfloat << std::setfill('0') << std::setw(2) << N
           << "." << std::setfill('0') << std::setw(2) << Lt << "x"
           << std::setfill('0') << std::setw(2) << Ls << ".beta" << std::setw(5)
           << std::fixed << std::setprecision(3) << beta << "_eps"
           << epsstr.str() << "." << std::setfill('0') << std::setw(7)
           << udint(rd) % (10000000) << start_cond << integr_str << ".dat";

      fname = outs.str();
      outs.str("");

      outs << "# N=" << N << "\n"
           << "# Lt=" << Lt << "\n"
           << "# Ls=" << Ls << "\n"
           << "# beta=" << beta << "\n"
           << "# eps=" << eps_t << "\n"
           << "# thermalisation sweeps=" << thermal << "\n"
           << "# measurement    sweeps=" << swps << "\n#\n"
           << "# Measurment\tEnergyDensity\tMeanPlaq\txi_G\tchi_m\tQ\n";


      ofile.open(fname, std::ios::trunc);

      ofile << outs.str();
      ofile.close();
      outs.str("");

      outs << measurements;
      auto aux = outs.str();
      m_length = static_cast<int>(aux.length());
      outs.str("");
    }


    for (size_t i = 0; i < measurements; ++i) {

      if (!quiet) {
        std::cout << "\t\t(" << std::setw(m_length) << (i + 1) << "/"
                  << std::setw(m_length) << measurements << ")\t";
      }

      // Do Sweeps
      start = clck::now();
      for (auto ii = 0LU; ii < swps; ++ii) {
        if (integrator > 0) {
          rkint.update(Force, lattice, rlx);
        }
        else {
          euler.update(Force, lattice, rlx);
        }
      }
      end = clck::now();
      time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count();
      ;

      if (not isCPConfig(lattice, 1.e-10)) {
        std::cout
          << "WARNING: Configuration does not fulfill CP(N-1) constraints!!";
        auto errors = maxCPConfigError(lattice);
        std::cout << "Max Errors: (" << std::scientific << errors[0] << ","
                  << errors[1] << ")\n";
      }
      if (not no_meas) {
        auto energy = cp_act.energyDensity(cpn_bc);
        auto plaq = cp_obs.meanPlaquette(u1_bc);
        auto xi = cp_obs.xi_G();
        auto chi = cp_obs.chi_m();
        auto Q = cp_obs.Q_top(cpn_bc);
        outs << std::setfill(' ') << std::setprecision(6) << std::setw(10)
             << std::fixed << (i + 1) * swps << std::scientific << "\t"
             << std::setw(10) << energy << "\t" << std::setw(13) << plaq << "\t"
             << std::setw(13) << xi << "\t" << std::setw(13) << chi << "\t"
             << std::setprecision(2) << std::setw(6) << std::fixed << Q
             << std::endl;

        ofile.open(fname, std::ios::app);
        ofile << outs.str();
        ofile.close();

        if (!quiet) {
          std::cout << outs.str();
        }
        outs.str("");

        v_energy.push_back(energy);
        v_plaquette.push_back(plaq);
        v_xi_G.push_back(xi);
        v_chi_m.push_back(chi);
        v_Q.push_back(Q);
      }


      if (not outfile.empty()) {
        Hdf5File h5file(outfile);

        std::stringstream grp_str, lat_str;
        grp_str << "beta" << std::setfill('0') << std::setw(6) << std::fixed
                << std::setprecision(3) << beta;

        h5file.push(grp_str.str());

        lat_str << "cfg_" << std::setfill('0') << std::setw(4) << (i + 1);

        h5file.write_lattice(lattice, lat_str.str());

        // Not implemented yet
        // h5file.flush();

        // Write some measured observables as attributes
        if (not no_meas) {
          h5file.addAttribute("energy_density", v_energy.back(), lat_str.str());
          h5file.addAttribute("mean_plaquette", v_plaquette.back(),
                              lat_str.str());
        }
      }
    }
    auto finish = clck::now();


    double cfg_per_sec = 1.e9 * (static_cast<double>(thermal) +
                                 static_cast<double>(swps * measurements));
    cfg_per_sec /= (static_cast<double>(time));

    double total_time =
      static_cast<double>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(finish - begin)
          .count()) *
      1.e-9;


    auto mean_energy = vec_mean(v_energy);
    auto std_err_energy = vec_std_err(v_energy);
    auto mean_plaq = vec_mean(v_plaquette);
    auto std_err_plaq = vec_std_err(v_plaquette);
    auto mean_Q = vec_mean(v_Q);
    auto std_err_Q = vec_std_err(v_Q);
    auto chi_t = vec_exp2(v_Q) / static_cast<double>(lattice.volume());
    auto std_err_chi_t =
      vec_exp2_std_err(v_Q) / static_cast<double>(lattice.volume());
    auto mean_xi = vec_mean(v_xi_G);
    auto std_err_xi = vec_std_err(v_xi_G);
    auto mean_chi = vec_mean(v_chi_m);
    auto std_err_chi = vec_std_err(v_chi_m);


    std::cout << "\t\tEnergy density: " << mean_energy << "\t" << std_err_energy
              << "\n";
    std::cout << "\t\tPlaquette:      " << mean_plaq << "\t" << std_err_plaq
              << "\n";
    std::cout << "\t\tCorr. Length:   " << mean_xi << "\t" << std_err_xi
              << "\n";
    std::cout << "\t\tMag. Suscept.:  " << mean_chi << "\t" << std_err_chi
              << "\n";
    std::cout << "\t\tTop. charge:  " << mean_Q << "\t" << std_err_Q << "\n";
    std::cout << "\t\tTop. Suscept.:  " << chi_t << "\t" << std_err_chi_t
              << std::endl;

    // Write res to file
    resfile.open(resstr, std::ios::app);

    // beta mean_energy std_error_energy mean_plaq std_err_plaq ...
    // std_err_chit_t
    outs << std::setfill(' ') << std::fixed << std::setprecision(4)
         << std::setw(8) << beta << "\t" << std::scientific
         << std::setprecision(6) << std::setw(13) << mean_energy << "\t"
         << std::setw(13) << std_err_energy << "\t" << std::setw(13)
         << mean_plaq << "\t" << std::setw(13) << std_err_plaq << "\t"
         << std::setw(13) << mean_xi << "\t" << std::setw(13) << std_err_xi
         << "\t" << std::setw(13) << mean_chi << "\t" << std::setw(13)
         << std_err_chi << "\t" << std::setw(13) << mean_Q << "\t"
         << std::setw(13) << std_err_Q << "\t" << std::setw(13) << chi_t << "\t"
         << std::setw(13) << std_err_chi_t << "\t" << std::endl;

    resfile << outs.str();
    resfile.close();
    outs.str("");

    std::cout << "\tLangevin run for beta=" << beta << " finished."
              << std::endl;
    std::cout << "\tThis took " << total_time << " sec." << std::endl;
    std::cout << "\tWe generated " << cfg_per_sec << " configurations per sec."
              << std::endl;
    std::cout << "\t..........................................................."
                 "..........."
              << std::endl;
  }


  return EXIT_SUCCESS;
}

// ######################################################################
// Helper Functions
// ######################################################################

/// Initialise CP(N-1) Lattice with random sites and links
template <typename Generator>
void canonical_CP_random_init(CPlattice &lattice, Generator &gen)
{


  for (auto &&latit :
       lattice) //= lattice.begin(); latit != lattice.end(); ++latit)
  {

    // Random site
    latit.site = randomCP<N>(gen);

    // Random links
    for (auto i = 0lu; i < lattice.dimensions(); ++i) {
      latit.links[i] = randomU1(gen);
    }
  }
}


/// Check if constraints are fulfilled
bool isCPConfig(CPlattice const &lattice, double const &eps)
{

  bool res = true;

  for (auto l : lattice) {

    res = std::abs(l.site.norm() - 1) < eps;
    // std::cout << l.site << std::endl;
    if (not res) {
      return res;
    }
    for (auto s : l.links) {
      res = std::abs(std::abs(s.value()) - 1) < eps;
      if (not res) {
        return res;
      }
    }
  }


  return res;
}

std::vector<double> maxCPConfigError(CPlattice const &lattice)
{

  std::vector<double> res;

  res.push_back(0.);
  res.push_back(0.);

  for (auto l : lattice) {

    auto new_cp = std::abs(l.site.norm() - 1);
    if (new_cp > res[0]) {
      res[0] = new_cp;
    }

    for (auto s : l.links) {
      auto new_gauge = std::abs(std::abs(s.value()) - 1);
      if (new_gauge > res[1]) {
        res[1] = new_gauge;
      }
    }
  }

  return res;
}

#else  // CBLAS available
int main(int argc, char *argv[]) { return 0; }
#endif // CBLAS available
