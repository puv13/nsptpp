///! \file
/// A canonical MC implementation for the CP(N-1) model
///

#include "./canonical_cpn_mc.h"
#include "io.h"
#include <omp.h>

// Global constants known at compile time
const std::size_t dim = 2;

int main(int argc, char *argv[])
{

  // ----------------------------------------------------------------------
  // Say Hi
  // ----------------------------------------------------------------------

  std::cout << "\nRunning MC simulation of the CP(N-1) model ..." << std::endl;
  std::cout << "________________________________________________"
            << "______________________________\n"
            << std::endl;

  // ----------------------------------------------------------------------
  // Get MC parameters from command line options
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
  // Sanity check
  // ----------------------------------------------------------------------
  if (outfile.empty() && no_meas) {
    std::cout << "*************************************************************"
                 "*********\n";
    std::cout << "You don't want to do measurements and "
              << "you didn't specify a file\nto safe the "
              << "configurations.\n"
              << "Why even bother?\n"
              << "... I am out ... " << std::endl;
    std::cout << "*************************************************************"
                 "*********\n";

    return 1;
  }

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
  // std::size_t udint_max=udint.max();
  std::vector<size_t> seeds;
  // for (auto i=0ul; i<std::mt19937_64::state_size; ++i){
  // seeds.push_back(udint(rd));}
  // DEBUG
  for (auto i = 0ul; i < std::mt19937_64::state_size; ++i) {
    seeds.push_back(i);
  }

  std::seed_seq sseq(seeds.begin(), seeds.end());
  std::mt19937_64 generator(sseq);


  // ----------------------------------------------------------------------
  // Init CP(N)
  // ----------------------------------------------------------------------

  std::array<size_t, dim> dimsar;
  dimsar.fill(Ls);
  dimsar[0] = Lt;
  CP<std::complex<double>, N> cp_init(1.0, true);
  U1 u1_init(std::complex<double>(1.0, 1.0), true);

  // Boundary conditions
  // "Gauge field"
  std::array<double, dim> u1_bc_ar;
  u1_bc_ar.fill(1.0);
  BoundaryCondition<double, dim> u1_bc(u1_bc_ar, "PERIODIC U1");
  std::string cpn_str = "PERIODIC CPN";


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
  if (0 != twisted) {
    start_cond = start_cond + "_twist";
  }
  std::ofstream resfile;
  std::stringstream outs;
  outs << "N" << std::setfill('0') << std::setw(2) << N << "."
       << std::setfill('0') << std::setw(2) << Lt << "x" << std::setfill('0')
       << std::setw(2) << Ls << start_cond << ".dat";
  auto resstr = outs.str();
  outs.str("");
  resfile.open(resstr, std::ios::app);
  resfile << "\n\n# "
             "beta\tEnergyDensity\tError\tMeanPlaq\tError\txi_G\tError\tchi_"
             "m\tError\tQ\tError\tchi_t\tError\tAccCPN\tAccU1\n";
  resfile.close();

  // ----------------------------------------------------------------------
  // Loop over beta values
  // ----------------------------------------------------------------------
  std::cout << "Looping over beta ..." << std::endl;
  std::cout << "_______________________________________________________________"
               "_______________\n"
            << std::endl;


  for (auto beta : scan_range) {
    // ----------------------------------------------------------------------
    // MC
    // ----------------------------------------------------------------------
    std::cout << "\n\tWorking on beta=" << beta << std::endl;
    std::cout << "\t..........................................................."
                 "..........."
              << std::endl;
    auto begin = clck::now();

    // Update epsilons
    double cp_eps = (beta < 1) ? 3.5 : 1. / beta;
    double u1_eps = (beta < 1) ? 0.95 : 0.5;

    // Update, Action and Observable init
    FullLattice<CP<std::complex<double>, N>, U1, dim> cplat(dimsar, cp_init,
                                                            u1_init);

    canonical_CP_update<std::complex<double>, N, dim, decltype(no_twist)>
      MC_update(&cplat, beta, cpn_bc, u1_bc);
    CanonicalCPNAction<std::complex<double>, N, dim> cp_act(cplat, beta);
    CanonicalCPNObservables<std::complex<double>, N, dim> cp_obs(cplat, beta);

    // Hot or cold start?
    if (hot_start) {
      canonical_CP_random_init(cplat, generator);
    }

    // Tune acceptance ?
    if (tune_accept) {

      // Tuning parameters
      double acc_diff = 0.025;
      double acc_goal = 0.6;
      size_t acc_sweep = 20;
      size_t acc_max_iter = 100;

      std::cout << "\tTuning acceptance rates ... " << std::endl;
      // "Thermalise" to make sure acceptance rates are stable
      auto acceptance =
        MC_update.Sweep(generator, acc_sweep * 10, cp_eps, u1_eps);

      auto ca = acceptance[0];
      auto ua = acceptance[1];

      auto count = 0lu;

      // Init parameters
      bool cp_conv = std::abs(ca - acc_goal) < acc_diff;
      bool u1_conv = std::abs(ua - acc_goal) < acc_diff;
      double cp_fact = 1.;
      double u1_fact = 1.;

      while ((not(cp_conv and u1_conv)) and (count < acc_max_iter)) {
        // Don't go on for ever if suitable parameters cannot be found
        count++;

        // Modify CP update
        if (not cp_conv) {
          (ca - acc_goal) < 0 ? cp_fact *= 0.9 : cp_fact *= 1.1;
        }

        // Modify U1 update
        if (not u1_conv) {
          (ua - acc_goal) < 0 ? u1_fact *= 0.9 : u1_fact *= 1.1;
        }

        acceptance = MC_update.Sweep(generator, 20 * N, cp_eps * cp_fact,
                                     u1_eps * u1_fact);

        ca = acceptance[0];
        ua = acceptance[1];

        cp_conv = std::abs(ca - acc_goal) < acc_diff;
        u1_conv = std::abs(ua - acc_goal) < acc_diff;

        if (!quiet) {
          std::cout << std::fixed << std::setprecision(3);
          std::cout << "\t\t"
                    << "iteration: " << std::setw(6) << count << "\t"
                    << "cp_eps: " << std::setw(8) << cp_eps * cp_fact << "\t"
                    << "cp_acc: " << std::setw(4) << ca << "\t"
                    << "u1_eps: " << std::setw(8) << u1_eps * u1_fact << "\t"
                    << "u1_acc: " << std::setw(4) << ua << std::endl;
        }
      }

      std::cout << "\n\tFinished tuning" << std::endl;
      if (std::abs(ca - acc_goal) > acc_diff) {
        std::cout << "\tWARNING: Could not find suitable parameters for cp_eps"
                  << std::endl;
      }
      if (std::abs(ua - acc_goal) > acc_diff) {
        std::cout << "\tWARNING: Could not find suitable parameters for u1_eps"
                  << std::endl;
      }

      cp_eps *= cp_fact;
      u1_eps *= u1_fact;

      std::cout << "\tUsing cp_eps = " << cp_eps << " and u1_eps = " << u1_eps
                << std::endl
                << std::endl;
    }


    // Start MC run
    std::cout << "\tStarting MC run with (Lt,Ls)=(" << Lt << "," << Ls << ")"
              << ", N=" << N << " and beta=" << beta << std::endl;
    // std::cout <<
    // "\t......................................................................"
    // << std::endl;

    // Thermalisation
    std::cout << "\tThermalisation ... " << std::endl;
    auto start = clck::now();

    // Show progress of thermalisation
    const size_t prog_steps = 100;
    progress_bar pb(static_cast<double>(prog_steps), 56, std::string("\t\t"));
    for (auto t = 0LU; t < prog_steps; ++t) {
      pb.progress(static_cast<double>(t));
      size_t th = static_cast<size_t>(std::ceil(
        static_cast<double>(thermal) / static_cast<double>(prog_steps)));
      MC_update.Sweep(generator, th, cp_eps, u1_eps);
    }
    pb.progress(static_cast<double>(prog_steps));
    std::cout << std::endl;


    auto end = clck::now();
    auto time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "\t\t... done\n" << std::endl;

    // Vectors to store results of measurements
    std::vector<double> v_energy;
    std::vector<double> v_plaquette;
    std::vector<double> v_xi_G;
    std::vector<double> v_chi_m;
    std::vector<double> v_acc_cpn;
    std::vector<double> v_acc_u1;
    std::vector<double> v_Q;

    int m_length;
    std::ofstream ofile;
    std::string fname;
    if (not no_meas) {

      // Starting Measurements
      std::cout << "\tMeasurements ... " << std::endl;


      // Filename and header for file output
      outs.str("");
      outs << "N" << std::defaultfloat << std::setfill('0') << std::setw(2) << N
           << "." << std::setfill('0') << std::setw(2) << Lt << "x"
           << std::setfill('0') << std::setw(2) << Ls << ".beta" << std::setw(5)
           << std::fixed << std::setprecision(3) << beta << "."
           << std::setfill('0') << std::setw(7) << udint(rd) % (10000000)
           << start_cond << ".dat";

      fname = outs.str();
      outs.str("");

      outs << "# N=" << N << "\n"
           << "# Lt=" << Lt << "\n"
           << "# Ls=" << Ls << "\n"
           << "# beta=" << beta << "\n"
           << "# cp_eps=" << cp_eps << "\n"
           << "# u1_eps=" << u1_eps << "\n"
           << "# thermalisation sweeps=" << thermal << "\n"
           << "# measurement    sweeps=" << swps << "\n#\n"
           << "# "
              "Measurment\tEnergyDensity\tMeanPlaq\txi_G\tchi_"
              "m\tQ\tAcceptCP\tAcceptU1\n";


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
      auto acceptance = MC_update.Sweep(generator, swps, cp_eps, u1_eps);
      end = clck::now();
      time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count();
      ;

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
             << std::setprecision(2) << std::setw(6) << std::fixed << Q << "\t"
             << "acc:" << std::scientific << std::setw(7) << acceptance[0]
             << "\t" << std::setprecision(2) << std::setw(7) << acceptance[1]
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
        v_acc_cpn.push_back(acceptance[0]);
        v_acc_u1.push_back(acceptance[1]);
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

        h5file.write_lattice(cplat, lat_str.str());

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
    auto mean_Q = vec_mean(v_Q);
    auto std_err_Q = vec_std_err(v_Q);
    auto chi_t = vec_exp2(v_Q) / static_cast<double>(cplat.volume());
    auto std_err_chi_t =
      vec_exp2_std_err(v_Q) / static_cast<double>(cplat.volume());
    auto mean_plaq = vec_mean(v_plaquette);
    auto std_err_plaq = vec_std_err(v_plaquette);
    auto mean_xi = vec_mean(v_xi_G);
    auto std_err_xi = vec_std_err(v_xi_G);
    auto mean_chi = vec_mean(v_chi_m);
    auto std_err_chi = vec_std_err(v_chi_m);

    auto bs_xi = bootstrap(v_xi_G, 10000, generator);

    std::cout << "\t\tEnergy density: " << mean_energy << "\t" << std_err_energy
              << "\n";
    std::cout << "\t\tPlaquette:      " << mean_plaq << "\t" << std_err_plaq
              << "\n";
    std::cout << "\t\tCorr. Length:   " << mean_xi << "\t" << std_err_xi
              << "\n";
    std::cout << "\t\tCorr. Length BS:" << bs_xi[0] << "\t" << bs_xi[1] << "\n";
    std::cout << "\t\tMag. Suscept.:  " << mean_chi << "\t" << std_err_chi
              << "\n";
    std::cout << "\t\tTop. charge:  " << mean_Q << "\t" << std_err_Q << "\n";
    std::cout << "\t\tTop. Suscept.:  " << chi_t << "\t" << std_err_chi_t
              << std::endl;

    // for (auto i=2LU; i<=2048; i*=2)
    // {

    //     std::cout << "Bin size: " << i << "\tMean" <<
    //     vec_mean(binning(v_xi_G,i))
    // 	      << "\tError: " << vec_std_err(binning(v_xi_G,i)) << std::endl;
    // }

    // Write res to file
    resfile.open(resstr, std::ios::app);

    // beta mean_energy std_error_energy mean_plaq std_err_plaq ...
    // mean_acc_cpn mean_acc_u1
    outs << std::setfill(' ') << std::fixed << std::setprecision(4)
         << std::setw(8) << beta << "\t" << std::scientific
         << std::setprecision(6) << std::setw(13) << mean_energy << "\t"
         << std::setw(13) << std_err_energy << "\t" << std::setw(13)
         << mean_plaq << "\t" << std::setw(13) << std_err_plaq << "\t"
         << std::setw(13) << mean_xi << "\t" << std::setw(13) << std_err_xi
         << "\t" << std::setw(13) << mean_chi << "\t" << std::setw(13)
         << std_err_chi << "\t" << std::setw(13) << mean_Q << "\t"
         << std::setw(13) << std_err_Q << "\t" << std::setw(13) << chi_t << "\t"
         << std::setw(13) << std_err_chi_t << "\t" << std::setw(13)
         << vec_mean(v_acc_cpn) << "\t" << std::setw(13) << vec_mean(v_acc_u1)
         << std::endl;

    resfile << outs.str();
    resfile.close();
    outs.str("");

    std::cout << "\tMC run for beta=" << beta << " finished." << std::endl;
    std::cout << "\tThis took " << total_time << " sec." << std::endl;
    std::cout << "\tWe generated " << cfg_per_sec << " configurations per sec."
              << std::endl;
    std::cout << "\t..........................................................."
                 "..........."
              << std::endl;
  }
  std::cout << "\n... Done" << std::endl;
  std::cout << "_______________________________________________________________"
               "_______________\n"
            << std::endl;
  return 0;
}


// ######################################################################
// Helper Functions
// ######################################################################

/// Initialise CP(N-1) Lattice with random sites and links
template <typename FullLatticeType, typename generator>
void canonical_CP_random_init(FullLatticeType &lattice, generator &gen)
{

  for (auto &&latit :
       lattice) //= lattice.begin(); latit != lattice.end(); ++latit)
  {
    // Random site
    latit.site = randomCP<latit.site.size()>(gen);

    // Random links
    for (auto i = 0lu; i < lattice.dimensions(); ++i) {
      latit.links[i] = randomU1(gen);
    }
  }
}
