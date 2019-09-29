///! \file
/// Header for the implementation of NSPT for the CP(N-1) model
///

#pragma once

#include "../src/actions/canonical_cpn.h"
#include "../src/expansion.h"
#include "../src/gaugegroups/sun.h"
#include "../src/gaugegroups/u1.h"
#include "../src/integrator.h"
#include "../src/integrators/euler.h"
#include "../src/integrators/rk.h"
#include "../src/lattice.h"
#include "../src/latticefields/cp.h"
#include "../src/progress_bar.h"
#include "../src/stat.h"
#include "cpn_nspt_params.h"
#include "io.h"
#include "randgen.h"
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
const std::size_t dim = 2;

#ifdef SET_ORDER_
constexpr size_t order = SET_ORDER_;
#else
constexpr size_t order = 7;
#endif

constexpr size_t TH = 100;

// Using
// ----------------------------------------------------------------------
using namespace CanonicalCPNExpHelpers;
using namespace CommonCPNExpHelpers;

using cplx = std::complex<double>;
using CPN = CP<cplx, N>;

using CPlattice = FullLattice<Expansion<CPN, order>, Expansion<U1, order>, dim>;

using clck = std::chrono::high_resolution_clock;

template <size_t O>
using CPL = FullLattice<Expansion<CPN, O>, Expansion<U1, O>, dim>;


// Declarations
//
template <typename Generator>
void Lattice_Init(CPlattice &lattice, Generator &gen);

template <size_t ord, typename P = double>
auto EnergyDensity(CPL<ord> const &lat, const BoundaryCondition<P, dim> &sbc =
                                          BoundaryCondition<P, dim>())
  -> Expansion<cplx, ord>;


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
        auto cor = -1. * (nrm[i] - aux) / aux;

        site[i] *= cor;
        // norm has changed !!!
        nrm = site * site;
      }
      else {
        std::cout << "WARNING: CAN NOT NORMALISE CPN" << std::endl;
      }
    }
  }
}

template <size_t ord>
void normalise(Expansion<U1, ord> &link, double unit_eps = 1.e-10)
{

  // auto A = log(link);

  // A=real(A);
  // A=(A-dagger(A))*0.5;


  // link = exp(A);

  auto unit = dagger(link) * link;

  // First order is special
  if (std::abs(unit[0].value() - 1.0) > unit_eps) {
    link[0] *= std::sqrt(unit[0].value());
  }

  // Higher orders
  for (auto i = 1LU; i < ord; ++i) {
    if (std::abs(unit[i].value()) > unit_eps) {

      auto aux = 2. * std::real((dagger(link[0]) * link[i]).value());
      if (0.0 != aux) {
        auto cor = -1. * (unit[i].value() - aux) / aux;

        link[i] *= cor;
        // absolute value has changed !!!
        unit = dagger(link) * link;
      }
      else {
        std::cout << "WARNING: CAN NOT NORMALISE U1" << std::endl;
      }
    }
  }
}

template <size_t ord>
void NormaliseLattice(
  FullLattice<Expansion<CPN, ord>, Expansion<U1, ord>, dim> &lattice,
  double eps = 1.e-14)
{

  for (auto &&latit : lattice) {
    // Site
    normalise<ord>(latit.site, eps);

    // Links
    for (auto mu = 0lu; mu < lattice.dimensions(); ++mu) {
      normalise<ord>(latit.links[mu], eps);
    }
  }
}


/// Bring phase to (-pi, pi] interval
/// How to do this efficiently is not entirely trivial.
/// See, e.g., https://stackoverflow.com/a/11980755 
inline U1 wrapPhase(U1 const &phase)
{

  double PI2 = M_PI * 2.0;
  double val = std::imag(phase.value());

  val = val - PI2 * floor(val / PI2);

  val = (val > M_PI) ? val - M_PI : val;


  return U1(val) * cplx(0., 1.);
}


template <size_t ord>
void StochasticGaugeFixing(CPL<ord> &lattice, double const &gauge_eps)
{
  std::vector<Expansion<U1, ord>> gauge_phases;
  gauge_phases.resize(lattice.volume());


  // First Sweep
  // Compute the phases of the gauge transformation
  for (auto latit = lattice.begin(); latit != lattice.end(); ++latit) {
    Expansion<U1, ord> phase;

    auto U_x = (*latit).links;

    for (auto mu = 0LU; mu < lattice.dimensions(); ++mu) {

      auto U_x_mu = latit.neighborLink(mu, mu, Direction::FORWARD);
      // Sign is intentional!!
      phase = phase + (log(U_x[mu]) - log(U_x_mu));
    }
    gauge_phases[latit.index()] = phase;

    // Debug:
    // std::cout << "idx: "  << latit.index();
    for (auto o = 0LU; o < ord; ++o) {
      gauge_phases[latit.index()][o] =
        wrapPhase(gauge_phases[latit.index()][o]);
      // std::cout << "\t" << gauge_phases[latit.index()][o] << std::endl;
    }
  }


  // Second Sweep
  // Apply gauge transformation
  for (auto latit = lattice.begin(); latit != lattice.end(); ++latit) {

    auto left = exp(gauge_eps * gauge_phases[latit.index()]);
    (*latit).site = (*latit).site * left;

    for (auto mu = 0LU; mu < lattice.dimensions(); ++mu) {
      auto right =
        exp(-1.0 * gauge_eps * gauge_phases[latit.neighbor(mu).index()]);

      (*latit).links[mu] = left * (*latit).links[mu] * right;
    }
  }
}

template <size_t ord>
void ZeroModeSubtraction(CPL<ord> &lattice)
{


  std::array<Expansion<U1, ord>, dim> ZeroMode;
  // First Sweep
  // Compute the zero mode contribution
  for (auto latit = lattice.begin(); latit != lattice.end(); ++latit) {
    for (auto mu = 0LU; mu < lattice.dimensions(); ++mu) {
      ZeroMode[mu] = ZeroMode[mu] + log((*latit).links[mu]);
    }
  }

  for (auto mu = 0LU; mu < lattice.dimensions(); ++mu) {
    ZeroMode[mu] = ZeroMode[mu] * (1. / static_cast<double>(lattice.volume()));
  }

  // Second Sweep
  // Remove the zero mode contribution
  for (auto latit = lattice.begin(); latit != lattice.end(); ++latit) {
    for (auto mu = 0LU; mu < lattice.dimensions(); ++mu) {

      auto A = log((*latit).links[mu]);

      A = A - ZeroMode[mu];

      // While we are at it ...
      for (auto o = 0LU; o < ord; ++o) {
        A[o] = wrapPhase(A[o]);
      }
      // A[0]=U1(0.0);

      (*latit).links[mu] = exp(A);
    }
  }
}


/// Thermalise a lattice order by order in the expansion
template <size_t ord>
CPL<ord> Sequential_Thermalisation(ranlux24_normal<double> &randgen,
                                   std::ostream &outs, int int_type = 1)
{
  std::array<size_t, dim> dimsar;
  dimsar.fill(Ls);
  dimsar[0] = Lt;
  CPN cp_init(1.0, true);
  U1 u1_init(1.0);
  CPL<ord> res(dimsar, cp_init, u1_init);

  // Recursion
  // ----------------------------------------------------------------------
  auto lower = Sequential_Thermalisation<ord - 1>(randgen, outs, int_type);
  for (auto i = 0LU; i < res.volume(); ++i) {
    res[i].site = lower[i].site.template cast_higher<ord>();
    // Normalise up to new order
    // ----------------------------------------------------------------------
    auto corr = res[i].site[1] * res[i].site[ord - 1 - 1];
    for (auto l = 2LU; l < ord - 1; l++) {
      corr += res[i].site[l] * res[i].site[ord - 1 - l];
    }

    auto t = 0LU;
    for (t = 0LU; t < N; ++t) {
      if (0. != std::real(res[i].site[0][t])) {
        corr = -1. * corr / (2. * std::real(res[i].site[0][t]));
        res[i].site[ord - 1][t] = corr;
        break;
      }
      throw std::runtime_error("Cannot normalise field.");
    }

    // Normalise up to new order
    // ----------------------------------------------------------------------
    for (auto mu = 0lu; mu < res.dimensions(); ++mu) {
      res[i].links[mu] = lower[i].links[mu].template cast_higher<ord>();

      auto corrl = dagger(res[i].links[mu][1]) * res[i].links[mu][ord - 1 - 1];

      for (auto l = 2LU; l < ord - 1; l++) {
        corrl += dagger(res[i].links[mu][l]) * res[i].links[mu][ord - 1 - l];
      }
      corrl = corrl * (-1. / (2. * std::real(res[i].links[mu][0].value())));

      res[i].links[mu][ord - 1] = corrl;
    }
  }

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


  // Thermalisation
  // ----------------------------------------------------------------------
  Expan::CPForce<N, dim, ord, CPtwistPhase<std::complex<double>, N>, double>
    myforce(cpn_bc, u1_bc);
  IntegratorsNoExp::RKInt<
    CPL<ord>,
    std::vector<CanonicalCPNExpHelpers::CPForceStructExp<N, dim, ord>>,
    CPCasimirStruct<N>, decltype(myforce),
    randomCPdf<N, dim, ranlux24_normal<double>>, ranlux24_normal<double>>
    rkint;
  rkint.setBeta(1.);
  rkint.setStepSize(eps_t);

  IntegratorsNoExp::EulerInt<
    CPL<ord>,
    std::vector<CanonicalCPNExpHelpers::CPForceStructExp<N, dim, ord>>,
    decltype(myforce), randomCPdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>
    euler;
  euler.setBeta(1.);
  euler.setStepSize(eps_t);

  outs << "\n\n";
  for (auto i = 0LU; i < TH; ++i) {

    // Log progress
    auto ED = EnergyDensity<ord>(res);
    // auto dev=OrderByOrderDeviationFromConstraints<ord>(res);

    outs << ord << "\t" << i << std::setfill(' ') << std::setprecision(6);
    for (auto k = 0LU; k < ord; ++k) {
      outs << std::scientific << "\t" << std::setw(16) << std::real(ED[k]);
    }

    // for (auto o=0LU; o<ord; ++o){
    //    outs << "\t"  << std::setw(16) << dev[o];
    // }
    outs << std::endl;


    if (int_type > 0) {
      rkint.update(myforce, res, randgen);
    }
    else {
      euler.update(myforce, res, randgen);
    }

    // StochasticGaugeFixing<ord>(res,eps_t);
    // ZeroModeSubtraction<ord>(res);
  }

  // NormaliseLattice<ord>(res);


  return res;
}


/// Thermalisation recursion stops at order 1
const size_t first_ord = 1LU;
template <>
CPL<first_ord>
Sequential_Thermalisation<first_ord>(ranlux24_normal<double> &randgen,
                                     std::ostream &outs, int int_type)
{
  const size_t ord = first_ord;
  std::array<size_t, dim> dimsar;
  dimsar.fill(Ls);
  dimsar[0] = Lt;
  CPN cp_init(1.0, true);
  U1 u1_init(1.0);
  CPL<ord> res(dimsar, cp_init, u1_init);

  // Set lowest order
  // ----------------------------------------------------------------------
  for (auto &&latit : res) {
    // Site
    latit.site[0] = CP<cplx, N>(1., true);

    // Links
    for (auto i = 0lu; i < res.dimensions(); ++i) {
      latit.links[i][0] = U1(1.0);
    }

    // StochasticGaugeFixing<first_ord>(res,eps_t);
    // ZeroModeSubtraction<first_ord>(res);
  }


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


  // Thermalisation
  // ----------------------------------------------------------------------
  Expan::CPForce<N, dim, ord, CPtwistPhase<std::complex<double>, N>, double>
    myforce(cpn_bc, u1_bc);
  IntegratorsNoExp::RKInt<
    CPL<ord>,
    std::vector<CanonicalCPNExpHelpers::CPForceStructExp<N, dim, ord>>,
    CPCasimirStruct<N>, decltype(myforce),
    randomCPdf<N, dim, ranlux24_normal<double>>, ranlux24_normal<double>>
    rkint;
  rkint.setBeta(1.);
  rkint.setStepSize(eps_t);

  IntegratorsNoExp::EulerInt<
    CPL<ord>,
    std::vector<CanonicalCPNExpHelpers::CPForceStructExp<N, dim, ord>>,
    decltype(myforce), randomCPdf<N, dim, ranlux24_normal<double>>,
    ranlux24_normal<double>>
    euler;
  euler.setBeta(1.);
  euler.setStepSize(eps_t);


  outs << "\n\n";
  for (auto i = 0LU; i < TH; ++i) {
    // Log progress
    auto ED = EnergyDensity<ord>(res);
    // auto dev=OrderByOrderDeviationFromConstraints<ord>(res);

    outs << ord << "\t" << i << std::setfill(' ') << std::setprecision(6);
    for (auto k = 0LU; k < ord; ++k) {
      outs << std::scientific << "\t" << std::setw(16) << std::real(ED[k]);
    }

    // for (auto o=0LU; o<ord; ++o){
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

  // NormaliseLattice<ord>(res);

  // StochasticGaugeFixing<ord>(res,eps_t);
  // ZeroModeSubtraction<ord>(res);


  return res;
}


std::array<double, 4> DeviationFromConstraints(CPlattice const &lattice)
{

  std::array<double, 4> dev;
  dev.fill(0.0);
  RunningStatistics<double> rss;
  RunningStatistics<double> rsl;

  for (const auto &latit : lattice) {
    // Site
    auto site = latit.site;
    auto sitenrm = site * site;


    double deviation = std::abs(sitenrm[0] - 1.0);

    for (auto i = 1LU; i < order; ++i) {
      deviation += std::abs(sitenrm[i]);
    }


    rss.push(deviation);
    dev[0] = rss.mean();
    dev[1] = (dev[1] < deviation) ? deviation : dev[1];

    // Links
    for (auto mu = 0lu; mu < lattice.dimensions(); ++mu) {
      auto linknrm = dagger(latit.links[mu]) * latit.links[mu];

      double linkdev = std::abs(linknrm[0].value() - 1.0);

      for (auto i = 1LU; i < order; ++i) {
        linkdev += std::abs(linknrm[i].value());
      }

      rsl.push(linkdev);
      dev[3] = (dev[3] < linkdev) ? linkdev : dev[3];
    }
    dev[2] = rsl.mean();
  }
  return dev;
}

template <size_t ord>
std::array<double, 2 * ord>
OrderByOrderDeviationFromConstraints(CPL<ord> const &lattice)
{

  std::array<double, 2 * ord> dev;
  dev.fill(0.0);
  RunningStatisticsArray<double, ord> rss;
  RunningStatisticsArray<double, ord> rsl;

  for (const auto &latit : lattice) {
    // Site
    auto site = latit.site;
    auto sitenrm = site * site;

    std::valarray<double> svs(0., ord);
    std::valarray<double> svl(0., ord);

    double deviation = std::abs(sitenrm[0] - 1.0);

    svs[0] = deviation;


    for (auto i = 1LU; i < ord; ++i) {
      svs[i] = std::abs(sitenrm[i]);
      deviation += svs[i];
    }

    rss.push(svs);


    // Links
    for (auto mu = 0lu; mu < lattice.dimensions(); ++mu) {
      auto linknrm = dagger(latit.links[mu]) * latit.links[mu];

      double linkdev = std::abs(linknrm[0].value() - 1.0);
      svl[0] = linkdev;

      for (auto i = 1LU; i < ord; ++i) {
        svl[i] = std::abs(linknrm[i].value());
        linkdev += svl[i];
      }

      rsl.push(svl);
      // dev[3]= (dev[3] < linkdev) ? linkdev : dev[3];
    }
    // dev[2]=rsl.mean();
  }

  for (auto i = 0LU; i < ord; ++i) {
    dev[i] = rss.mean()[i];
    dev[ord + i] = rsl.mean()[i];
  }


  return dev;
}
