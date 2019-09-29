///! \file
/// Header for the implementation of NSPT for the CP(N-1) model with quartic
/// action
///

#pragma once

#include "actions/quartic_cpn.h"
#include "expansion.h"
#include "gaugegroups/sun.h"
#include "integrator.h"
#include "integrators/euler.h"
#include "integrators/rk.h"
#include "io.h"
#include "lattice.h"
#include "latticefields/cp.h"
#include "progress_bar.h"
#include "quartic_cpn_nspt_params.h"
#include "randgen.h"
#include "stat.h"
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
// ----------------------------------------------------------------------
constexpr std::size_t dim = 2;
constexpr std::size_t TH = 100;


#ifdef SET_ORDER_
constexpr size_t order = SET_ORDER_;
#else
constexpr size_t order = 7;
#endif


// Using
// ----------------------------------------------------------------------
using namespace QuarticCPNExpHelpers;
using namespace QuarticCPNHelpers;
using namespace CommonCPNExpHelpers;

using CPlattice = SiteLattice<Expansion<CPN, order>, dim>;
using clck = std::chrono::high_resolution_clock;

template <size_t O>
using CPL = SiteLattice<Expansion<CPN, O>, dim>;

template <size_t ord, typename P = double>
auto QuarticEnergyDensity(
  CPL<ord> const &lat,
  const BoundaryCondition<P, dim> &bc = BoundaryCondition<P, dim>())
  -> Expansion<cplx, ord>;

template <typename Generator>
void Lattice_Init(CPlattice &lattice, Generator &gen);


template <size_t ord>
std::array<double, ord>
OrderByOrderDeviationFromConstraints(CPL<ord> const &lattice)
{

  std::array<double, ord> dev;
  dev.fill(0.0);
  RunningStatisticsArray<double, ord> rss;

  for (const auto &latit : lattice) {
    auto site = latit;
    auto sitenrm = site * site;

    std::valarray<double> svs(0., ord);

    double deviation = std::abs(sitenrm[0] - 1.0);

    svs[0] = deviation;

    for (auto i = 1LU; i < ord; ++i) {
      svs[i] = std::abs(sitenrm[i]);
      deviation += svs[i];
    }

    rss.push(svs);
  }

  for (auto i = 0LU; i < ord; ++i) {
    dev[i] = rss.mean()[i];
  }


  return dev;
}


template <size_t ord>
void normalise(Expansion<CPN, ord> &site, double nrm_eps = 1.e-10)
{

  auto nrm = site * site;

  // First order is special
  if (std::abs(nrm[0] - 1.0) > nrm_eps) {
    site[0] *= std::sqrt(nrm[0]);
  }

  // Higher orders
  for (auto i = 1LU; i < ord; ++i) {
    if (std::abs(nrm[i]) > nrm_eps) {
      auto aux = 2. * std::real(site[0] * site[i]);
      if (0.0 != aux) {

        double cor_alt_den = std::real(nrm[i] - aux);
        double sgn = cor_alt_den < 0. ? -1.0 : 1.0;
        double sgn_aux = aux < 0. ? -1.0 : 1.0;

        auto cor = -1. * (nrm[i] - aux) / aux;

        cor_alt_den *= sgn;
        aux *= sgn_aux;

        // auto aux_alt =
        // -1.*sgn*sgn_aux*std::exp(std::log(cor_alt_den)-std::log(aux));


        // std::cout << "COR: "  << std::abs(1.0-std::real(cor)) << "\t" <<
        // std::abs(1.0-aux_alt)
        // 	  << std::endl;

        site[i] *= std::real(cor);
        // norm has changed !!!
        nrm = site * site;
      }
      else {
        std::cout << "WARNING: CAN NOT NORMALISE CPN" << std::endl;
      }
    }
  }
}


/// Thermalise a lattice order by order in the expansion
template <size_t ord>
CPL<ord> Sequential_Thermalisation(ranlux24_normal<double> &randgen,
                                   int int_type = 1)
{
  std::array<size_t, dim> dimsar;
  dimsar.fill(Ls);
  dimsar[0] = Lt;
  CPL<ord> res(dimsar);

  // Recursion
  // ----------------------------------------------------------------------
  auto lower = Sequential_Thermalisation<ord - 1>(randgen, int_type);
  for (auto i = 0LU; i < res.volume(); ++i) {
    res[i] = lower[i].template cast_higher<ord>();
    // Normalise up to new order
    // ----------------------------------------------------------------------
    auto corr = res[i][1] * res[i][ord - 1 - 1];
    for (auto l = 2LU; l < ord - 1; l++) {
      corr += res[i][l] * res[i][ord - 1 - l];
    }

    auto t = 0LU;
    for (t = 0LU; t < N; ++t) {
      if (0. != std::real(res[i][0][t])) {
        corr = -1. * corr / (2. * std::real(res[i][0][t]));
        res[i][ord - 1][t] = corr;
        break;
      }
      throw std::runtime_error("Cannot normalise field.");
    }
  }

  // ----------------------------------------------------------------------
  // Init Boundaries
  // ----------------------------------------------------------------------

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

  std::string cpn_str = "PERIODIC CPN";
  BoundaryCondition<CPtwistPhase<std::complex<double>, N>, dim> cpn_bc(tp_ar,
                                                                       cpn_str);


  // Thermalisation
  // ----------------------------------------------------------------------
  Expan::QuarticCPForce<N, dim, ord, CPtwistPhase<std::complex<double>, N>>
    myforce(cpn_bc);

  IntegratorsNoExp::RKInt<
    CPL<ord>, typename std::vector<Expansion<nummat<N>, ord>>,
    QuarticCPCasimirStruct<N>, decltype(myforce),
    QuarticCPNExpHelpers::QuarticrandomCPdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>
    rkint;
  rkint.setBeta(1.);
  rkint.setStepSize(eps_t);

  IntegratorsNoExp::EulerInt<
    CPL<ord>, typename std::vector<Expansion<nummat<N>, ord>>,
    decltype(myforce),
    QuarticCPNExpHelpers::QuarticrandomCPdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>
    euler;
  euler.setBeta(1.);
  euler.setStepSize(eps_t);


  for (auto i = 0LU; i < TH; ++i) {
    if (int_type > 0) {
      rkint.update(myforce, res, randgen);
    }
    else {
      euler.update(myforce, res, randgen);
    }
  }


  // Normalise again
  const double nrm_eps = 1.e-14;


  for (auto i = 0LU; i < res.volume(); ++i) {
    // Normalise highest order
    // ----------------------------------------------------------------------


    auto nrm = res[i] * res[i];

    // std::cout << std::endl;
    // for (auto j=0LU; j<ord; ++j)
    // {
    //     std::cout << std::abs(nrm[j]) << "\t";
    // }
    // std::cout << std::endl;

    normalise(res[i], nrm_eps);


    nrm = res[i] * res[i];

    // std::cout << std::endl;
    // for (auto j=0LU; j<ord; ++j)
    // {
    //     std::cout << std::abs(nrm[j]) << "\t";
    // }
    // std::cout << std::endl << std::endl;


    // if (std::abs(nrm[ord-1]) > nrm_eps){
    //     auto aux = 2.*std::real(res[i][0]*res[i][ord-1]);
    //     if (0.0 != aux ){
    // 	auto cor = -1.*(nrm[ord-1]-aux)/aux;

    // 	res[i][ord-1]*=cor;
    // 	// norm has changed !!!
    // 	nrm=res[i]*res[i];

    // 	for (auto j=0LU; j<ord; ++j)
    // 	{
    // 	    std::cout << std::abs(nrm[j]) << "\t";
    // 	}
    // 	std::cout << std::endl << std::endl;

    //     }
    //     else{
    // 	std::cout << "WARNING: CAN NOT NORMALISE CPN" << std::endl;
    //     }
    // }
  }

  return res;
}


/// Thermalisation recursion stops at order 1
const size_t first_ord = 1LU;
template <>
CPL<first_ord>
Sequential_Thermalisation<first_ord>(ranlux24_normal<double> &randgen,
                                     int int_type)
{
  const size_t ord = first_ord;
  std::array<size_t, dim> dimsar;
  dimsar.fill(Ls);
  dimsar[0] = Lt;
  CPN cp_init(1.0, true);
  CPL<ord> res(dimsar, cp_init);

  // Set lowest order
  // ----------------------------------------------------------------------
  for (auto &&latit : res) {
    latit[0] = CP<cplx, N>(1., true);
  }


  // ----------------------------------------------------------------------
  // Init Boundaries
  // ----------------------------------------------------------------------

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

  std::string cpn_str = "PERIODIC CPN";
  BoundaryCondition<CPtwistPhase<std::complex<double>, N>, dim> cpn_bc(tp_ar,
                                                                       cpn_str);


  // Thermalisation
  // ----------------------------------------------------------------------
  Expan::QuarticCPForce<N, dim, ord, CPtwistPhase<std::complex<double>, N>>
    myforce(cpn_bc);

  IntegratorsNoExp::RKInt<
    CPL<ord>, typename std::vector<Expansion<nummat<N>, ord>>,
    QuarticCPCasimirStruct<N>, decltype(myforce),
    QuarticCPNExpHelpers::QuarticrandomCPdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>
    rkint;
  rkint.setBeta(1.);
  rkint.setStepSize(eps_t);

  IntegratorsNoExp::EulerInt<
    CPL<ord>, typename std::vector<Expansion<nummat<N>, ord>>,
    decltype(myforce),
    QuarticCPNExpHelpers::QuarticrandomCPdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>
    euler;
  euler.setBeta(1.);
  euler.setStepSize(eps_t);


  for (auto i = 0LU; i < TH; ++i) {
    if (int_type > 0) {
      rkint.update(myforce, res, randgen);
    }
    else {
      euler.update(myforce, res, randgen);
    }
  }

  return res;
}
