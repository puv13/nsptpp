///! \file
/// A MC implementation for compact QED
///

#include "./u1_mc.h"
#include "io.h"
#include <omp.h>

// Global constants known at compile time
#ifdef SET_DIM_
const std::size_t dim = SET_DIM_;
#else
const std::size_t dim = 4;
#endif


int main(int argc, char *argv[])
{
  // ----------------------------------------------------------------------
  // Say Hi
  // ----------------------------------------------------------------------
  std::cout << "\nRunning MC simulation of compact QED ..." << std::endl;
  std::cout << "________________________________________________"
            << "______________________________\n"
            << std::endl;

  // ----------------------------------------------------------------------
  // Get MC parameters from command line options
  // ----------------------------------------------------------------------

  parse_command_line_options(argc, argv);
  print_params();

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
  // Default beta range
  // ----------------------------------------------------------------------
  std::vector<double> scan_range;
  if (beta0 <= 0.0) {

    auto num_beta = 8;
    for (int i = 0; i <= num_beta; ++i) {
      auto high = 1.2;
      auto low = 0.8;

      auto delta = high - low;

      scan_range.push_back(low + delta * static_cast<double>(i) / num_beta);
    }
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
  std::cout << std::setw(28) << "\tdim:" << dim << std::endl;
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
  for (auto i = 0ul; i < std::mt19937_64::state_size; ++i) {
    seeds.push_back(udint(rd));
  }

  std::seed_seq sseq(seeds.begin(), seeds.end());
  std::mt19937_64 generator(sseq);


  // ----------------------------------------------------------------------
  // Init Lattice
  // ----------------------------------------------------------------------

  // Lattice dimensions
  std::array<size_t, dim> dimsar;
  dimsar.fill(Ls);
  dimsar[0] = Lt;

  // U1
  U1 u1_init(std::complex<double>(1.0, 1.0), true);

  // Boundary conditions
  std::array<double, dim> u1_bc_ar;
  u1_bc_ar.fill(0.0);
  BoundaryCondition<double, dim> u1_bc(u1_bc_ar, "PERIODIC U1");

  // ----------------------------------------------------------------------
  // Resfile header
  // ----------------------------------------------------------------------
  std::string start_cond = (hot_start) ? "" : "_cs";

  std::ofstream resfile;
  std::stringstream outs;
  outs << "cQED" << dim << "." << std::setfill('0') << std::setw(2) << Lt << "x"
       << std::setfill('0') << std::setw(2) << Ls << start_cond << ".dat";
  auto resstr = outs.str();
  outs.str("");
  resfile.open(resstr, std::ios::app);
  resfile << "\n\n# beta\tEnergyDensity\tError\tAccU1\n";
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

    // Update epsilon
    double u1_eps = 1;

    // Lattice and Action init
    LinkLattice<U1, dim> qplat(dimsar, u1_init);
    QEDAction<dim> QED_act(qplat, beta);
    Compact_QED_update<dim> mc(&qplat, beta);

    // Hot or cold start?
    if (hot_start) {
      mc.prepare_hot(generator);
    }

    // std::cout.precision(6);
    // std::cout << std::fixed  <<  QED_act.energyDensity() << std::endl;


    // Tune acceptance ?
    // overrelaxation is switched off for acceptance tuning
    if (tune_accept) {

      // Tuning parameters
      double acc_diff = 0.025;
      double acc_goal = 0.6;
      size_t acc_sweep = 20;
      size_t acc_max_iter = 100;

      // Sweep(generator &gen, const double &eps =0.5, size_t nsweep=1, size_t
      // nhit=1, size_t nor=1)
      std::cout << "\tTuning acceptance rates ... " << std::endl;

      // or_param = 0
      auto acceptance =
        mc.Sweep(generator, u1_eps, acc_sweep * 10, hit_param, 0);

      auto count = 0lu;

      bool u1_conv = std::abs(acceptance - acc_goal) < acc_diff;
      double u1_fact = 1.;

      while ((not u1_conv) and (count < acc_max_iter)) {
        // Don't go on for ever if suitable parameters cannot be found
        count++;

        // Modify U1 update
        (acceptance - acc_goal) < 0 ? u1_fact *= 0.9 : u1_fact *= 1.1;

        acceptance = mc.Sweep(generator, u1_eps * u1_fact, 20, hit_param, 0);
        u1_conv = std::abs(acceptance - acc_goal) < acc_diff;

        if (!quiet) {
          std::cout << std::fixed << std::setprecision(3);
          std::cout << "\t\t"
                    << "iteration: " << std::setw(6) << count << "\t"
                    << "u1_eps: " << std::setw(8) << u1_eps * u1_fact << "\t"
                    << "u1_acc: " << std::setw(4) << acceptance << std::endl;
        }
      }

      std::cout << "\n\tFinished tuning" << std::endl;
      if (std::abs(acceptance - acc_goal) > acc_diff) {
        std::cout << "\tWARNING: Could not find suitable parameters for u1_eps"
                  << std::endl;
      }

      u1_eps *= u1_fact;

      std::cout << "\tUsing  u1_eps = " << u1_eps << std::endl << std::endl;
    }


    // Start MC run
    std::cout << "\tStarting MC run with  beta=" << beta << std::endl;
    // std::cout <<
    // "\t......................................................................"
    // << std::endl;

    // Thermalisation
    // Overrelaxation is switched off for thermalisation
    std::cout << "\tThermalisation ... " << std::endl;
    auto start = clck::now();
    // Sweep(generator &gen, const double &eps =0.5, size_t nsweep=1, size_t
    // nhit=1, size_t nor=1)

    // Show progress of thermalisation
    const size_t prog_steps = 100;
    progress_bar pb(static_cast<double>(prog_steps), 56, std::string("\t\t"));
    for (auto t = 0LU; t < prog_steps; ++t) {
      pb.progress(static_cast<double>(t));
      size_t th = static_cast<size_t>(std::ceil(
        static_cast<double>(thermal) / static_cast<double>(prog_steps)));

      mc.Sweep(generator, u1_eps, th, hit_param, 0);
    }
    pb.progress(static_cast<double>(prog_steps));
    std::cout << std::endl;


    auto end = clck::now();
    auto time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "\t\t... done\n" << std::endl;

    // Vectors to store results of measurements
    std::vector<double> v_energy;
    std::vector<double> v_acc_u1;

    int m_length;
    std::ofstream ofile;
    std::string fname;

    if (not no_meas) {

      // Starting Measurements
      std::cout << "\tMeasurements ... " << std::endl;


      // Filename and header for file output
      outs.str("");
      outs << "cQED" << std::defaultfloat << dim << "." << std::setfill('0')
           << std::setw(2) << Lt << "x" << std::setfill('0') << std::setw(2)
           << Ls << ".beta" << std::setw(5) << std::fixed
           << std::setprecision(3) << beta << "." << std::setfill('0')
           << std::setw(7) << udint(rd) % (10000000) << start_cond << ".dat";

      fname = outs.str();
      outs.str("");

      outs << "# dim=" << dim << "\n"
           << "# Lt=" << Lt << "\n"
           << "# Ls=" << Ls << "\n"
           << "# beta=" << beta << "\n"
           << "# u1_eps=" << u1_eps << "\n"
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
      // Sweep(generator &gen, const double &eps =0.5, size_t nsweep=1, size_t
      // nhit=1, size_t nor=1)
      auto acceptance = mc.Sweep(generator, u1_eps, swps, hit_param, or_param);
      end = clck::now();
      time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count();
      ;

      if (not no_meas) {
        auto energy = QED_act.energyDensity(u1_bc);
        outs << std::setfill(' ') << std::setprecision(6) << std::setw(10)
             << std::fixed << (i + 1) * swps << std::scientific << "\t"
             << std::setw(10) << energy << "\t"
             << "acc:" << std::scientific << std::setw(7) << acceptance
             << std::endl;

        ofile.open(fname, std::ios::app);
        ofile << outs.str();
        ofile.close();

        if (!quiet) {
          std::cout << outs.str();
        }
        outs.str("");

        v_energy.push_back(energy);
        v_acc_u1.push_back(acceptance);
      }


      if (not outfile.empty()) {
        Hdf5File h5file(outfile);

        std::stringstream grp_str, lat_str;
        grp_str << "beta" << std::setfill('0') << std::setw(6) << std::fixed
                << std::setprecision(3) << beta;

        h5file.push(grp_str.str());

        lat_str << "cfg_" << std::setfill('0') << std::setw(4) << (i + 1);

        h5file.write_lattice(qplat, lat_str.str());

        // Not implemented yet
        // h5file.flush();

        // Write some measured observables as attributes
        if (not no_meas) {
          h5file.addAttribute("energy_density", v_energy.back(), lat_str.str());
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

    std::cout << "\t\tEnergy density: " << mean_energy << "\t" << std_err_energy
              << "\n";

    // Write res to file
    resfile.open(resstr, std::ios::app);

    // beta mean_energy std_error_energy mean_plaq std_err_plaq ...
    // mean_acc_cpn mean_acc_u1
    outs << std::setfill(' ') << std::fixed << std::setprecision(4)
         << std::setw(8) << beta << "\t" << std::scientific
         << std::setprecision(6) << std::setw(13) << mean_energy << "\t"
         << std::setw(13) << std_err_energy << "\t" << std::setw(13)
         << vec_mean(v_acc_u1) << std::endl;

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
