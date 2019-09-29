///! \file
/// Header for the implementation of the Langevin algorithm for the CP(N-1)
/// model
/// with quartic action

#pragma once

#include "../src/actions/quartic_cpn.h"
#include "../src/gaugegroups/sun.h"
#include "../src/gaugegroups/u1.h"
#include "../src/integrator.h"
#include "../src/integrators/euler.h"
#include "../src/integrators/rk.h"
#include "../src/lattice.h"
#include "../src/latticefields/cp.h"
#include "../src/progress_bar.h"
#include "../src/stat.h"
#include "quartic_cpn_langevin_params.h"
#include <algorithm>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>


// Global constants known at compile time
constexpr std::size_t dim = 2;


using namespace QuarticCPNHelpers;

using clck = std::chrono::high_resolution_clock;
using QuarticCPlattice = SiteLattice<CPN, dim>;

template <typename Generator>
void quartic_CP_random_init(QuarticCPlattice &lattice, Generator &gen);

bool isQuarticCPConfig(QuarticCPlattice const &lattice,
                       double const &eps = 1.e-10);
double maxQuarticCPConfigError(QuarticCPlattice const &lattice);

template <typename Generator>
class QuarticrandomCPdf {
 public:
  nummat<N> operator()(Generator &rng) const
  {
    nummat<N> rd;
    rd = noiseSUN<N, Generator>(rng, CP_Langevin_noise_scale);
    return rd;
  }
};
