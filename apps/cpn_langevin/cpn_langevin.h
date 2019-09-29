///! \file
/// Header for the implementation of the Langevin algorithm for the CP(N-1)
/// model
///

#pragma once

#include "../src/actions/canonical_cpn.h"
#include "../src/gaugegroups/sun.h"
#include "../src/gaugegroups/u1.h"
#include "../src/integrator.h"
#include "../src/integrators/euler.h"
#include "../src/integrators/rk.h"
#include "../src/lattice.h"
#include "../src/latticefields/cp.h"
#include "../src/progress_bar.h"
#include "../src/stat.h"
#include "cpn_langevin_params.h"
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

using namespace CanonicalCPNHelpers;

using clck = std::chrono::high_resolution_clock;
using CPlattice = FullLattice<CPN, U1, dim>;

template <typename Generator>
void canonical_CP_random_init(CPlattice &lattice, Generator &gen);

bool isCPConfig(CPlattice const &lattice, double const &eps = 1.e-10);
std::vector<double> maxCPConfigError(CPlattice const &lattice);


template <typename Generator>
class randomCPdf {
 public:
  CPForceStruct<N, dim> operator()(Generator &rng) const
  {

    CPForceStruct<N, dim> rd;
    rd.site = noiseSUN<N, Generator>(rng, CP_Langevin_noise_scale);


    for (auto i = 0LU; i < rd.links.size(); ++i) {
      rd.links[i] = noiseU1(rng, CP_Langevin_noise_scale);
    }

    return rd;
  }
};
