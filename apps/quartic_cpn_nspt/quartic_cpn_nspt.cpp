/// \file
///
/// Numerical stochastic perturbation theory for the CP(N-1) model with quartic
/// action

#include "quartic_cpn_nspt.h"

int main(int argc, char *argv[])
{

  // ----------------------------------------------------------------------
  // Say Hi
  // ----------------------------------------------------------------------

  std::cout << "\nRunning NSPT for the CP(N-1) model ..." << std::endl;
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
  std::cout << std::setw(28) << "\torder:" << order << std::endl;
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
  // Resfile header
  // ----------------------------------------------------------------------
  std::string start_cond = (seq_init) ? "" : "_no_seqinit";
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
       << std::setw(2) << Ls << "QUARTIC_ord" << std::setfill('0')
       << std::setw(2) << order << "_eps" << epsstr.str() << "." << of_id.str()
       << start_cond << integr_str << ".dat";
  auto resstr = outs.str();
  outs.str("");
  resfile.open(resstr, std::ios::app);
  resfile << "\n\n#  Meas.\tEnergy density coefficients\t Average dev. from "
             "norm\t Max dev from norm\n";
  resfile.close();


  // ----------------------------------------------------------------------
  // Init Boundaries
  // ----------------------------------------------------------------------

  // CPN
  CPtwistPhase<std::complex<double>, N> no_twist(std::complex<double>(1.0));
  std::array<CPtwistPhase<std::complex<double>, N>, dim> tp_ar;
  tp_ar.fill(no_twist);

  std::string cpn_str = "PERIODIC CPN";
  if (twisted == 1) {

    CPtwistPhase<std::complex<double>, N> twist(std::complex<double>(0.0));

    auto ii = std::complex<double>(0., 1.);
    for (auto i = 0LU; i < N; ++i) {
      twist[i] = std::exp(-ii * (std::atan(1.) * 8. / static_cast<double>(N)) *
                          static_cast<double>(i));
    }


    tp_ar[1] = twist;
    cpn_str = "TWISTED FIRST SPATIAL";
  }


  BoundaryCondition<CPtwistPhase<std::complex<double>, N>, dim> cpn_bc(tp_ar,
                                                                       cpn_str);


  // ----------------------------------------------------------------------
  // Init stuff
  // ----------------------------------------------------------------------
  ranlux24_normal<double> rlx;

  std::array<size_t, dim> dimsar;
  dimsar.fill(Ls);
  dimsar[0] = Lt;
  CPN cp_init(1.0, true);
  CPlattice lattice(dimsar, cp_init);
  Expan::QuarticCPForce<N, dim, order, CPtwistPhase<std::complex<double>, N>>
    myforce(cpn_bc);

  if (debug) {
    rlx.set_generator_state("1857670 9122615 9091361 7165633 6622926 5093178 "
                            "10658136 14316459 2246149 476568 820330 11401107 "
                            "5667087 14569279 13287542 9551677 2350107 409682 "
                            "3168448 8334688 10338889 15771630 12488960 "
                            "15078927 0 0 0");
    rlx.set_distribution_state(
      "0.00000000000000000e+00 1.00000000000000000e+00 0");
  }

  Lattice_Init(lattice, generator);


  using namespace QuarticCPNHelpers;
  IntegratorsNoExp::RKInt<
    CPlattice, typename std::vector<Expansion<nummat<N>, order>>,
    QuarticCPCasimirStruct<N>, decltype(myforce),
    QuarticCPNExpHelpers::QuarticrandomCPdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>
    rkint;
  rkint.setBeta(1.);
  rkint.setStepSize(eps_t);

  IntegratorsNoExp::EulerInt<
    CPlattice, typename std::vector<Expansion<nummat<N>, order>>,
    decltype(myforce),
    QuarticCPNExpHelpers::QuarticrandomCPdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>
    euler;
  euler.setBeta(1.);
  euler.setStepSize(eps_t);


  // Sequential Thermalisation
  // ----------------------------------------------------------------------

  if (seq_init && seed_file.empty()) {
    lattice = Sequential_Thermalisation<order>(rlx, integrator);
  }
  if (not seed_file.empty()) {
    Hdf5File h5file(seed_file, true);
    h5file.read_lattice(lattice, "seed");


    std::string bcstring;
    h5file.readAttribute(bcstring, "boundary_condition", "seed");

    if (use_seed_state) {
      std::string seedstate;
      h5file.readAttribute(seedstate, "state", "seed");
      rlx.set_generator_state(seedstate);
    }

    if (bcstring != cpn_bc.getName()) {
      throw std::runtime_error("Seed Lattice has the wrong boundary conditon!");
    }
  }


  // Thermalisation
  // ----------------------------------------------------------------------
  auto begin = clck::now();

  std::cout << "\tThermalisation ... " << std::endl;
  auto start = clck::now();

  // Show progress of thermalisation
  const size_t prog_steps = 100;
  progress_bar pb(static_cast<double>(prog_steps), 56, std::string("\t\t"));
  for (auto t = 0LU; t < prog_steps; ++t) {
    pb.progress(static_cast<double>(t));
    size_t th = static_cast<size_t>(std::ceil(static_cast<double>(thermal) /
                                              static_cast<double>(prog_steps)));
    for (auto i = 0LU; i < th; ++i) {
      if (integrator > 0) {
        rkint.update(myforce, lattice, rlx);
      }
      else {
        euler.update(myforce, lattice, rlx);
      }
    }
  }
  pb.progress(static_cast<double>(prog_steps));
  std::cout << std::endl;

  auto end = clck::now();
  auto time =
    std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
  std::cout << "\t\t... done\n" << std::endl;


  // ******************************************************************************************

  // Observables
  // Vectors to store results of measurements
  std::vector<std::vector<double>> energy_coeff;
  for (auto i = 0LU; i < order; ++i) {
    energy_coeff.push_back(std::vector<double>());
  }

  outs.str("");

  outs << measurements;
  auto aux = outs.str();
  int m_length = static_cast<int>(aux.length());
  outs.str("");

  for (auto i = 0LU; i < measurements; ++i) {

    if (!quiet) {
      std::cout << "\t\t(" << std::setw(m_length) << (i + 1) << "/"
                << std::setw(m_length) << measurements << ")\t";
    }

    // Do Sweeps
    start = clck::now();
    for (auto ii = 0LU; ii < swps; ++ii) {

      if (integrator > 0) {
        rkint.update(myforce, lattice, rlx);
      }
      else {
        euler.update(myforce, lattice, rlx);
      }

      // if (gaugefix){
      // 	StochasticGaugeFixing<order>(lattice,eps_t);
      // }
      // if (zeromodes){
      // 	ZeroModeSubtraction<order>(lattice);
      // }
    }
    end = clck::now();
    time +=
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    ;

    if (not no_meas) {

      auto ED = QuarticEnergyDensity<order>(lattice);
      for (auto k = 0LU; k < order; ++k) {
        energy_coeff[k].push_back(std::real(ED[k]));
      }

      outs << std::setfill(' ') << std::setprecision(6) << std::setw(10)
           << std::fixed << (i + 1) * swps;

      for (auto k = 0LU; k < order; ++k) {
        outs << std::scientific << "\t" << std::setw(16)
             << energy_coeff[k].back();
      }

      auto dev = OrderByOrderDeviationFromConstraints(lattice);

      for (auto o = 0LU; o < order; ++o) {
        outs << "\t" << dev[o];
      }

      outs << std::endl;

      resfile.open(resstr, std::ios::app);
      resfile << outs.str();
      resfile.close();


      if (!quiet) {
        std::cout << outs.str();
      }
      outs.str("");
    }

    if (not outfile.empty() and not seed_only) {
      Hdf5File h5file(outfile);

      std::stringstream grp_str, lat_str;
      grp_str << "eps" << std::setfill('0') << std::setw(12) << std::fixed
              << std::setprecision(10) << eps_t;

      h5file.push(grp_str.str());

      lat_str << "cfg_" << std::setfill('0') << std::setw(4) << (i + 1);

      h5file.write_lattice(lattice, lat_str.str());
    }

    // Last Lattice is seed
    if (not outfile.empty() and i == measurements - 1) {
      Hdf5File h5file(outfile);

      h5file.write_lattice(lattice, "seed");
      h5file.addAttribute<std::string>("state", rlx.get_generator_state(),
                                       "seed");
      h5file.addAttribute<std::string>("boundary_condition", cpn_bc.getName(),
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

  std::cout << "\n\n";
  for (auto k = 0LU; k < order; ++k) {
    outs << std::scientific << "\t"
         << "Energy_Coef[" << std::setw(2) << k << "]:\t" << std::setw(14)
         << vec_mean(energy_coeff[k]) << " +/- " << std::setw(14)
         << vec_std_err(energy_coeff[k]) << std::endl;
  }

  std::cout << outs.str();
  std::cout << std::setfill(' ') << std::fixed << std::setprecision(4);
  std::cout << "\tThis took " << total_time << " sec." << std::endl;
  std::cout << "\tWe generated " << cfg_per_sec << " configurations per sec."
            << std::endl;
  std::cout << "\t............................................................."
               "........."
            << std::endl;


  // ******************************************************************************************


  return EXIT_SUCCESS;
}

// ######################################################################
// Helper Functions
// ######################################################################

/// Initialise CP(N-1) Lattice with random sites
template <typename Generator>
void Lattice_Init(CPlattice &lattice, Generator &gen)
{

  for (auto &&latit :
       lattice) //= lattice.begin(); latit != lattice.end(); ++latit)
  {
    // Random site
    latit[1] = 0. * randomCP<N>(gen);
    latit[0] = CP<cplx, N>(1., true);
  }
}


template <size_t ord, typename P = double>
auto QuarticEnergyDensity(CPL<ord> const &lat,
                          const BoundaryCondition<P, dim> &bc)
  -> Expansion<cplx, ord>
{
  Expansion<cplx, ord> res;

  // Loop over all sites
  for (auto fli = lat.begin(); fli != lat.end(); ++fli) {
    auto site = fli.site();

    for (auto i = 0lu; i < dim; i++) {
      auto sp = scalar_prod(site, fli.neighborSite(i, Direction::FORWARD, bc));
      res += sp * conj(sp);
    }
  }

  return static_cast<double>(N) / static_cast<double>(lat.volume()) * res;
}
