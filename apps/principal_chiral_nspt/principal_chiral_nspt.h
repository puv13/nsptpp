/// \file
///
/// Header for the implementation of NSPT for the principal chiral model

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
#include "principal_chiral_nspt_params.h"
#include <algorithm>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <sstream>


// -----------------------------------------------------------------------------
// Global constants known at compile time
// -----------------------------------------------------------------------------
constexpr std::size_t dim = 2;
constexpr std::size_t TH = 500;


#ifdef SET_ORDER_
constexpr size_t order = SET_ORDER_;
#else
constexpr size_t order = 7;
#endif


using namespace PCMHelpers;
using namespace PCMExpHelpers;

using clck = std::chrono::high_resolution_clock;
using cplx = std::complex<double>;
using PCMlattice = SiteLattice<Expansion<nummat<N>, order>, dim>;


template <size_t O>
using PCML = SiteLattice<Expansion<nummat<N>, O>, dim>;

template <size_t O>
using PCML = SiteLattice<Expansion<nummat<N>, O>, dim>;


template <size_t ord, typename P = double>
auto PCMEnergyDensity(PCML<ord> const &lat,
                      const BoundaryCondition<P, dim> &bc =
                        BoundaryCondition<P, dim>()) -> Expansion<cplx, ord>;

/// Initialise Lattice
template <typename Generator>
void Lattice_Init(PCMlattice &lattice, Generator &gen);

/// Calculate the deviation from Unitarity
template <size_t ord>
std::array<double, ord>
OrderByOrderDeviationFromConstraints(PCML<ord> const &lattice)
{

  std::array<double, ord> dev;
  dev.fill(0.0);
  RunningStatisticsArray<double, ord> rss;

  for (const auto &latit : lattice) {
    auto site = latit;
    auto sitenrm = site * dagger(site);

    std::valarray<double> svs(0., ord);

    double deviation =
      std::abs(normF(sitenrm[0]) - std::sqrt(static_cast<double>(N))) /
      static_cast<double>(N * N);

    svs[0] = deviation;

    for (auto i = 1LU; i < ord; ++i) {
      svs[i] = normF(sitenrm[i]) / static_cast<double>(N * N);
      deviation += svs[i];
    }

    rss.push(svs);
  }

  for (auto i = 0LU; i < ord; ++i) {
    dev[i] = rss.mean()[i];
  }

  return dev;
}


/// Unitarise expansion
template <size_t ord>
void normalise(Expansion<nummat<N>, ord> &site)
{

  auto aux = log(site);

  for (auto i = 0LU; i < ord; ++i) {
    // Make aux Hermitean
    aux[i] = 0.5 * (aux[i] - dagger(aux[i]));
    // Make aux traceless
    aux[i] += nummat<N>(-1. / static_cast<double>(N)) * trace(aux[i]);
  }

  aux = exp(aux);

  site = aux;
}

/// Unitarise all lattice fields
template <size_t ord>
void normalise(PCML<ord> &lat)
{
  for (auto i = 0LU; i < lat.volume(); ++i)
  // for (auto &fli=lat.begin(); fli != lat.end(); ++fli)
  // for (auto fli : lat)
  {
    normalise(lat[i]);
  }
}


/// Order by order thermalisation of field coefficients
template <size_t ord>
void Sequential_Thermalisation(PCML<order> &res,
                               ranlux24_normal<double> &randgen,
                               std::ostream &outs, int int_type = 1)
{
  std::array<size_t, dim> dimsar;
  dimsar.fill(Ls);
  dimsar[0] = Lt;
  PCML<ord> *aux = new PCML<ord>(dimsar);

  // Cast down and normalise to given order
  for (auto i = 0LU; i < aux->volume(); ++i) {
    (*aux)[i] = res[i].template cast_lower<ord>();
    normalise((*aux)[i]);
  }


  // -------------------------------------------------------------------------
  // Init Boundaries
  // -------------------------------------------------------------------------

  PCMtwistPhase<std::complex<double>, N> no_twist(std::complex<double>(1.0));
  std::array<PCMtwistPhase<std::complex<double>, N>, dim> tp_ar;
  tp_ar.fill(no_twist);

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

    tp_ar[1] = twist;
  }

  std::string pcm_str = "PERIODIC PCM";
  BoundaryCondition<PCMtwistPhase<std::complex<double>, N>, dim> pcm_bc(
    tp_ar, pcm_str);


  // Thermalisation
  // -------------------------------------------------------------------------
  Expan::PCMForce<N, dim, ord, PCMtwistPhase<std::complex<double>, N>>
    *myforce =
      new Expan::PCMForce<N, dim, ord, PCMtwistPhase<std::complex<double>, N>>(
        pcm_bc);

  auto *rkint = new IntegratorsNoExp::RKInt<
    PCML<ord>, typename std::vector<Expansion<nummat<N>, ord>>,
    PCMCasimirStruct<N>, decltype(*myforce),
    PCMExpHelpers::randomPCMdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>;
  rkint->setBeta(1.);
  rkint->setStepSize(eps_t);

  auto *euler = new IntegratorsNoExp::EulerInt<
    PCML<ord>, typename std::vector<Expansion<nummat<N>, ord>>,
    decltype(*myforce),
    PCMExpHelpers::randomPCMdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>;
  euler->setBeta(1.);
  euler->setStepSize(eps_t);

  for (auto i = 0LU; i < TH; ++i) {
    // Log progress
    auto ED = PCMEnergyDensity<ord>(*aux, pcm_bc);
    // auto dev=OrderByOrderDeviationFromConstraints<ord>(*aux);

    outs << ord << "\t" << i << std::setfill(' ') << std::setprecision(6);
    for (auto k = 0LU; k < ord; ++k) {
      outs << std::scientific << "\t" << std::setw(16) << std::real(ED[k]);
    }

    // for (auto o=0LU; o<ord; ++o){
    //    outs << "\t"  << std::setw(16) << dev[o];
    // }
    outs << std::endl;

    if (int_type > 0) {
      rkint->update(*myforce, *aux, randgen);
    }
    else {
      euler->update(*myforce, *aux, randgen);
    }
  }
  delete myforce;
  delete rkint;
  delete euler;

  // Cast back up
  for (auto i = 0LU; i < res.volume(); ++i) {
    res[i] = (*aux)[i].template cast_higher<order>();
    // normalise((*aux)[i]);
  }

  delete aux;

  // std::cout << ord << std::endl;
  Sequential_Thermalisation<ord + 1>(res, randgen, outs, int_type);
}

/// Terminate recursion for sequential thermalisation
template <>
void Sequential_Thermalisation<order>(PCML<order> &res,
                                      ranlux24_normal<double> &randgen,
                                      std::ostream &outs, int int_type)
{

  normalise(res);
  // -------------------------------------------------------------------------
  // Init Boundaries
  // -------------------------------------------------------------------------

  PCMtwistPhase<std::complex<double>, N> no_twist(std::complex<double>(1.0));
  std::array<PCMtwistPhase<std::complex<double>, N>, dim> tp_ar;
  tp_ar.fill(no_twist);

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

    tp_ar[1] = twist;
  }

  std::string pcm_str = "PERIODIC PCM";
  BoundaryCondition<PCMtwistPhase<std::complex<double>, N>, dim> pcm_bc(
    tp_ar, pcm_str);

  // Thermalisation
  // -------------------------------------------------------------------------
  Expan::PCMForce<N, dim, order, PCMtwistPhase<std::complex<double>, N>>
    myforce(pcm_bc);

  IntegratorsNoExp::RKInt<
    PCML<order>, typename std::vector<Expansion<nummat<N>, order>>,
    PCMCasimirStruct<N>, decltype(myforce),
    PCMExpHelpers::randomPCMdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>
    rkint;
  rkint.setBeta(1.);
  rkint.setStepSize(eps_t);

  IntegratorsNoExp::EulerInt<
    PCML<order>, typename std::vector<Expansion<nummat<N>, order>>,
    decltype(myforce),
    PCMExpHelpers::randomPCMdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>
    euler;
  euler.setBeta(1.);
  euler.setStepSize(eps_t);

  for (auto i = 0LU; i < TH; ++i) {
    // Log progress
    auto ED = PCMEnergyDensity<order>(res, pcm_bc);
    // auto dev=OrderByOrderDeviationFromConstraints<order>(res);

    outs << order << "\t" << i << std::setfill(' ') << std::setprecision(6);
    for (auto k = 0LU; k < order; ++k) {
      outs << std::scientific << "\t" << std::setw(16) << std::real(ED[k]);
    }

    // for (auto o=0LU; o<order; ++o){
    //    outs << "\t"  << std::setw(16) << dev[o];
    // }
    outs << std::endl;

    if (int_type > 0) {
      rkint.update(myforce, res, randgen);
    }
    else {
      euler.update(myforce, res, randgen);
    }
  }
  normalise(res);
}
