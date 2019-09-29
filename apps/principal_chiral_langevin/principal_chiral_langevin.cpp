/// \file
///
/// Numerical stochastic perturbation theory for the CP(N-1) model

#include "principal_chiral_langevin.h"

int main(int argc, char *argv[])
{

  // ----------------------------------------------------------------------
  // Say Hi
  // ----------------------------------------------------------------------

  std::cout << "\nRunning Langevin for the principal chrial model ..."
            << std::endl;
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

  // ----------------------------------------------------------------------
  // Define betas by hand
  // ----------------------------------------------------------------------

  std::vector<double> scan_range;

  if (beta0 <= 0.0) {

    scan_range = std::vector<double>({ //.8, .85,.9,.95
                                       0.1, 0.125, 0.15, 0.175, 0.2, 0.225,
                                       0.25, 0.3, 0.4, 0.5, .8, .85, .9, .95 });
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

  PCMtwistPhase<std::complex<double>, N> no_twist(std::complex<double>(1.0));
  std::array<PCMtwistPhase<std::complex<double>, N>, dim> tp_ar;
  tp_ar.fill(no_twist);

  std::string pcm_str = "PERIODIC PCM";
  if (twisted == 1) {
    bool nu = (N % 2 == 0);

    auto ii = std::complex<double>(0., 1.);

    std::complex<double> factor = 1.0;

    if (nu) {
      factor = std::exp(ii * M_PI / static_cast<double>(N));
    }

    PCMtwistPhase<std::complex<double>, N> twist(factor);


    for (auto i = 0LU; i < N; ++i) {
      if (0 != i) {
        twist[i] *= std::exp(ii * 2. * M_PI * static_cast<double>(i) /
                             static_cast<double>(N));
      }
    }


    pcm_str = "TWISTED TEMPORAL";
    tp_ar[1] = twist;
  }

  BoundaryCondition<PCMtwistPhase<std::complex<double>, N>, dim> pcm_bc(
    tp_ar, pcm_str);


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
  std::stringstream of_id;
  file_identifier.empty()
    ? of_id << std::setfill('0') << std::setw(7) << udint(rd) % (10000000)
    : of_id << file_identifier;


  outs << "N" << std::setfill('0') << std::setw(2) << N << "."
       << std::setfill('0') << std::setw(2) << Lt << "x" << std::setfill('0')
       << std::setw(2) << Ls << "_PCM_eps" << epsstr.str() << "." << of_id.str()
       << start_cond << integr_str << ".dat";

  auto resstr = outs.str();
  outs.str("");
  resfile.open(resstr, std::ios::app);
  resfile << "\n\n# "
             "beta\tEnergyDensity\tError\tMeanPlaq\tError\txi_G\tError\tchi_"
             "m\tError\tQ\tError\tchi_t\tError\n";
  resfile.close();

  for (auto beta : scan_range) {

    std::cout << "\n\tWorking on beta=" << beta << std::endl;
    std::cout << "\t..........................................................."
                 "..........."
              << std::endl;
    auto begin = clck::now();

    // ----------------------------------------------------------------------
    // Langevin Driver
    // ----------------------------------------------------------------------
    nummat<N> pcm_init(1.0);

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
    PCMlattice lattice(dimsar, pcm_init);
    PCMAction<std::complex<double>, N, dim> pcm_act(lattice, beta);

    NoExpansion::PCMForce<N, dim, PCMtwistPhase<std::complex<double>, N>> Force(
      pcm_bc);


    IntegratorsNoExp::EulerInt<
      PCMlattice, typename std::vector<nummat<N>>, decltype(Force),
      PCMrandomCPdf<ranlux24_normal<double>>, ranlux24_normal<double>>
      euler;
    euler.setBeta(beta);
    euler.setStepSize(eps_t);

    IntegratorsNoExp::RKInt<PCMlattice, typename std::vector<nummat<N>>,
                            PCMCasimirStruct<N>, decltype(Force),
                            PCMrandomCPdf<ranlux24_normal<double>>,
                            ranlux24_normal<double>>
      rkint;
    rkint.setBeta(beta);
    rkint.setStepSize(eps_t);

    // Hot or cold start?
    if (hot_start) {
      PCM_random_init(lattice, generator);
    }


    if (not seed_file.empty()) {
      Hdf5File h5file(seed_file, true);

      std::stringstream grp_str, lat_str;

      grp_str << "eps" << std::setfill('0') << std::setw(12) << std::fixed
              << std::setprecision(10) << eps_t;

      h5file.push(grp_str.str());

      grp_str.str("");
      grp_str << "beta" << std::setfill('0') << std::setw(6) << std::fixed
              << std::setprecision(3) << beta;

      h5file.push(grp_str.str());


      std::string bcstring;
      h5file.readAttribute(bcstring, "boundary_condition", "seed");
      h5file.read_lattice(lattice, "seed");

      if (use_seed_state) {
        std::string seedstate;
        h5file.readAttribute(seedstate, "state", "seed");
        rlx.set_generator_state(seedstate);
      }

      if (bcstring != pcm_bc.getName()) {
        throw std::runtime_error(
          "Seed Lattice has the wrong boundary conditon!");
      }
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
    // 	std::vector<double> v_plaquette;
    // 	std::vector<double> v_xi_G;
    // 	std::vector<double> v_chi_m;
    // 	std::vector<double> v_Q;

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
           << std::fixed << std::setprecision(3) << beta << "_PCM_eps"
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
           << "# Measurment\tEnergyDensity\n";


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
      for (auto jj = 0LU; jj < swps; ++jj) {
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

      if (not isPCMConfig(lattice, 1.e-10)) {
        std::cout
          << "WARNING: Configuration does not fulfill CP(N-1) constraints!!";
        auto error = maxPCMConfigError(lattice);
        std::cout << "Max Error: " << std::scientific << error << "\n";
      }
      if (not no_meas) {
        auto energy = pcm_act.energyDensity(pcm_bc);
        outs << std::setfill(' ') << std::setprecision(6) << std::setw(10)
             << std::fixed << (i + 1) * swps << std::scientific << "\t"
             << std::setw(10) << energy << std::endl;

        ofile.open(fname, std::ios::app);
        ofile << outs.str();
        ofile.close();

        if (!quiet) {
          std::cout << outs.str();
        }
        outs.str("");

        v_energy.push_back(energy);
      }

      if (not outfile.empty() and not seed_only) {
        Hdf5File h5file(outfile);

        std::stringstream grp_str, lat_str;

        grp_str << "eps" << std::setfill('0') << std::setw(12) << std::fixed
                << std::setprecision(10) << eps_t;

        h5file.push(grp_str.str());

        grp_str.str("");
        grp_str << "beta" << std::setfill('0') << std::setw(6) << std::fixed
                << std::setprecision(3) << beta;

        h5file.push(grp_str.str());

        lat_str << "cfg_" << std::setfill('0') << std::setw(4) << (i + 1);

        h5file.write_lattice(lattice, lat_str.str());
      }

      // Last Lattice is seed
      if (not outfile.empty() and i == measurements - 1) {
        Hdf5File h5file(outfile);

        std::stringstream grp_str, lat_str;

        grp_str << "eps" << std::setfill('0') << std::setw(12) << std::fixed
                << std::setprecision(10) << eps_t;

        h5file.push(grp_str.str());

        grp_str.str("");
        grp_str << "beta" << std::setfill('0') << std::setw(6) << std::fixed
                << std::setprecision(3) << beta;

        h5file.push(grp_str.str());


        h5file.write_lattice(lattice, "seed");
        h5file.addAttribute<std::string>("state", rlx.get_generator_state(),
                                         "seed");
        h5file.addAttribute<std::string>("boundary_condition", pcm_bc.getName(),
                                         "seed");
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

    std::cout << std::scientific << std::setprecision(8) << std::setw(13)
              << "\t\tEnergy density: " << mean_energy << "\t" << std::setw(13)
              << std_err_energy << "\n";

    // Write res to file
    resfile.open(resstr, std::ios::app);

    // beta mean_energy std_error_energy mean_plaq std_err_plaq ...
    // std_err_chit_t
    outs << std::setfill(' ') << std::fixed << std::setprecision(6)
         << std::setw(8) << beta << "\t" << std::scientific
         << std::setprecision(8) << std::setw(13) << mean_energy << "\t"
         << std::setw(13) << std_err_energy << "\t" << std::endl;

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

  // ################################################################################


  return EXIT_SUCCESS;
}


// ######################################################################
// Helper Functions
// ######################################################################

/// Initialise PCM Lattice with random sites
template <typename Generator>
void PCM_random_init(PCMlattice &lattice, Generator &gen)
{
  for (auto &&latit :
       lattice) //= lattice.begin(); latit != lattice.end(); ++latit)
  {
    // Random site
    latit = randomSUN<N>(gen);
  }
}


// Check if constraints are fulfilled
bool isPCMConfig(PCMlattice const &lattice, double const &eps)
{

  bool res = true;
  for (auto l : lattice) {
    res = normF(l.mult_conj(l) - nummat<N>(1.)) < eps;
    // std::cout << l.site << std::endl;
    if (not res) {
      return res;
    }
  }


  return res;
}

double maxPCMConfigError(PCMlattice const &lattice)
{

  double res(0.0);
  for (auto l : lattice) {
    auto new_cp = normF(l.mult_conj(l) - nummat<N>(1.));
    if (new_cp > res) {
      res = new_cp;
    }
  }

  return res;
}
