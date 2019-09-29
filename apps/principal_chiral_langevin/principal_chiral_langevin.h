///! \file
/// Header for the implementation of Langevin for the principal chiral model
///

#pragma once

#include "../src/actions/principal_chiral_model.h"
#include "../src/expansion.h"
#include "../src/gaugegroups/sun.h"
#include "../src/integrator.h"
#include "../src/integrators/euler.h"
#include "../src/integrators/rk.h"
#include "../src/lattice.h"
#include "../src/progress_bar.h"
#include "../src/randgen.h"
#include "../src/stat.h"
#include "io.h"
#include "principal_chiral_langevin_params.h"
#include <algorithm>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <sstream>

// Global constants known at compile time
constexpr std::size_t dim = 2;


using namespace PCMHelpers;


using clck = std::chrono::high_resolution_clock;
using PCMlattice = SiteLattice<nummat<N>, dim>;

template <typename Generator>
void PCM_random_init(PCMlattice &lattice, Generator &gen);


template <typename Generator>
class PCMrandomCPdf {
 public:
  nummat<N> operator()(Generator &rng) const
  {
    nummat<N> rd;
    rd = noiseSUN<N, Generator>(rng, PCM_Langevin_noise_scale);
    return rd;
  }
};


bool isPCMConfig(PCMlattice const &lattice, double const &eps = 1.e-10);
double maxPCMConfigError(PCMlattice const &lattice);
