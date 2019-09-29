/// \file
///
/// \brief Naive compact QED action
/// \details The naice compact QED action implemented here does not know about
///  expansions.


#pragma once

#include "../gaugegroups/u1.h"
#include "../lattice.h"
#include "action.h"
#include "expansion.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <numeric>


#ifndef M_PI
#define M_PI 3.14159265358979323846 // Double precision pi
#endif


using size_t = std::size_t;


template <size_t dim>
class QEDAction {


 private:
  LinkLattice<U1, dim> const *lattice;
  double beta;


 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  // Default constructor
  QEDAction() : beta(1.0) { lattice = nullptr; }

  // Construct with given lattice
  QEDAction(LinkLattice<U1, dim> &orig, const double &coupling = 1.0)
    : lattice(&orig), beta(coupling)
  {
  }


  // **********************************************************************
  // Deconstructor
  // **********************************************************************

  ~QEDAction() = default;

  // **********************************************************************
  // Member Functions
  // **********************************************************************

  /// Setter for beta
  void setbeta(const double &coupling) { beta = coupling; }
  /// Getter for beta
  double getbeta() { return this->beta; }

  template <typename P = double>
  double atLink(
    size_t const &idx, size_t const &dir,
    const BoundaryCondition<P, dim> &bc = BoundaryCondition<P, dim>()) const
  {
    auto lat = *lattice;
    LinkLatticeIterator<const LinkLattice<U1, dim>> lli(idx, lat);

    if (idx >= lattice->volume()) {
      throw std::runtime_error("Site out of range");
    }
    if (dir >= lattice->dimensions()) {
      throw std::runtime_error("Direction out of range");
    }


    return beta * (static_cast<double>((dim - 1)) * 2. -
                   std::real((lli.link(dir)).value() * staples(idx, dir, bc)));
  }


  template <typename P = double>
  std::complex<double> staples(
    size_t const &idx, size_t const &dir,
    const BoundaryCondition<P, dim> &bc = BoundaryCondition<P, dim>()) const
  {


    auto lat = *lattice;
    LinkLatticeIterator<const LinkLattice<U1, dim>> lli(idx, lat);

    std::complex<double> A = std::complex<double>(0.0);
    for (auto nu = 0LU; nu < dim; ++nu) {
      if (nu != dir) {
        U1 plus(std::complex<double>(1., 0.));
        U1 minus(std::complex<double>(1., 0.));

        // U_ν(x+µ)dagger(U_µ(n+ν))dagger(U_ν(x))
        plus *= lli.neighborLink(dir, nu, Direction::FORWARD, bc);
        plus *= dagger(lli.neighborLink(nu, dir, Direction::FORWARD, bc));
        plus *= dagger(lli.link(nu));

        // access to site (x-ν)
        auto neib = lli.neighbor(nu, -1);

        // dagger(U_ν(x+µ-ν))dagger(U_µ(n-ν))U_ν(x-ν)
        minus *= dagger(neib.neighborLink(dir, nu, Direction::FORWARD, bc));
        minus *= dagger(neib.link(dir));
        minus *= neib.link(nu);

        A += plus.value() + minus.value();
      }
    }

    return A;
  }

  /// Returns (the real part of the trace of ) the plaquette
  template <typename P = double>
  auto plaquette(
    size_t site, size_t mu, size_t nu,
    BoundaryCondition<P, dim> const &bc = BoundaryCondition<P, dim>()) const
    -> decltype(std::real(trace(U1(1.))))
  {
    if (site >= lattice->volume()) {
      throw std::runtime_error("Site out of range");
    }
    if (mu >= dim) {
      throw std::runtime_error("mu out of range");
    }
    if (nu >= dim) {
      throw std::runtime_error("nu out of range");
    }
    if (mu == nu) {
      throw std::runtime_error("mu must not equal nu");
    }

    auto it = lattice->begin();
    it = it[site];

    // U_µν(x) =U_µ(x)U_ν(x+µ)dagger(U_µ(n+ν))dagger(U_ν(x))

    auto plaq = it.link(mu);
    plaq *= it.neighborLink(mu, nu, Direction::FORWARD, bc);
    plaq *= dagger(it.neighborLink(nu, mu, Direction::FORWARD, bc));
    plaq *= dagger(it.link(nu));

    return std::real(trace(plaq));
  }

  /// Compute mean plaquette in the mu-nu plane
  template <typename P = double>
  auto meanPlaquette(size_t mu, size_t nu, BoundaryCondition<P, dim> const &bc =
                                             BoundaryCondition<P, dim>()) const
    -> decltype(std::real(trace(U1(1.))))
  {

    auto res = plaquette(0lu, mu, nu, bc);
    for (auto i = 1ul; i < lattice->volume(); ++i) {
      res += plaquette(i, mu, nu);
    }

    return res / static_cast<double>(lattice->volume());
  }


  /// Compute global plaquette average
  template <typename P = double>
  auto meanPlaquette(
    BoundaryCondition<P, dim> const &bc = BoundaryCondition<P, dim>()) const
    -> decltype(std::real(trace(U1(1.))))
  {

    decltype(std::real(trace(U1(1.)))) res = 0.0;

    for (auto mu = 0lu; mu < dim; ++mu) {
      for (auto nu = mu + 1lu; nu < dim; ++nu) {
        auto res_mu_nu = plaquette(0ul, mu, nu);
        for (auto i = 1ul; i < lattice->volume(); ++i) {
          res_mu_nu += plaquette(i, mu, nu, bc);
        }

        res_mu_nu *= 1. / static_cast<double>(lattice->volume());

        res += res_mu_nu;
      }
    }

    return 2.0 * res / (static_cast<double>(dim * (dim - 1)));
  }


  template <typename P = double>
  auto energyDensity(
    BoundaryCondition<P, dim> const &bc = BoundaryCondition<P, dim>()) const
    -> decltype(std::real(trace(U1(1.))))
  {

    return (1. - meanPlaquette(bc));
  }


  auto GaugeFieldSpatialSum() -> std::vector<std::array<double, dim>>
  {

    auto dims = lattice->dimensionsArray();
    std::array<double, dim> A;
    std::vector<std::array<double, dim>> Avec;

    A.fill(0.);

    // Output vector has size of temporal lattice extend
    for (auto l = 0LU; l < dims[0]; ++l) {
      Avec.push_back(A); // Init with "zero"
    }


    // Loop over all lattice sites
    // std::vector<std::size_t> tvec;
    // tvec.resize(Avec.size());
    for (auto fli = lattice->begin(); fli != lattice->end(); ++fli) {
      auto cord = lattice->linearIndexToCoord(fli.index());
      auto t = cord[0];


      // tvec[t]++;
      for (auto mu = 0LU; mu < dim; ++mu) {
        (Avec[t])[mu] += std::sin((fli.link(mu)).phase());
      }
    }

    // No normalisation !!!
    return Avec;
  }
};

template <size_t dim, typename P = double>
class Compact_QED_update {

 private:
  LinkLattice<U1, dim> *const lattice;
  double beta;
  BoundaryCondition<P, dim> link_bc;


 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************
  Compact_QED_update()
    : lattice(nullptr), beta(1.0), link_bc(BoundaryCondition<double, dim>())
  {
  }

  Compact_QED_update(
    LinkLattice<U1, dim> *const l, double const &coupling,
    BoundaryCondition<P, dim> const &lb = BoundaryCondition<P, dim>())
    : lattice(l), beta(coupling), link_bc(lb)
  {
  }

  // **********************************************************************
  // Destructor
  // **********************************************************************

  ~Compact_QED_update() = default;


  // **********************************************************************
  // Member functions
  // **********************************************************************

  /// Helper function for MC update
  template <typename generator>
  bool accept(double const &res, generator &gen) const
  {
    // Always accept smaller new action
    if (res > 1.) {
      return true;
    }

    // MC for larger new action
    std::uniform_real_distribution<double> dist(0, 1.0);
    double rnd = dist(gen);

    return (rnd <= res);
  }

  /// Multi-hit MC
  template <typename generator>
  std::size_t multihit_MC(generator &gen, const std::size_t &hits = 1,
                          const double &eps = 1.e-2) const
  {

    std::size_t accepted = 0;

    // Iterate over lattice sites
    for (LinkLatticeIterator<LinkLattice<U1, dim>> it = this->lattice->begin();
         it != this->lattice->end(); ++it) {
      // Loop over links
      for (auto mu = 0LU; mu < this->lattice->dimensions(); ++mu) {
        // Get Link
        auto U_old = it.link(mu);

        // Compute staples
        std::complex<double> A = std::complex<double>(0.0);
        for (auto nu = 0LU; nu < this->lattice->dimensions(); ++nu) {
          if (nu != mu) {
            /// \todo Check boundary conditions in staples
            U1 plus(std::complex<double>(1., 0.));
            U1 minus(std::complex<double>(1., 0.));

            // U_ν(x+µ)dagger(U_µ(n+ν))dagger(U_ν(x))
            plus *= it.neighborLink(mu, nu, Direction::FORWARD, this->link_bc);
            plus *= dagger(
              it.neighborLink(nu, mu, Direction::FORWARD, this->link_bc));
            plus *= dagger(it.link(nu));

            // access to site (x-ν)
            auto neib = it.neighbor(nu, -1);

            // dagger(U_ν(x+µ-ν))dagger(U_µ(n-ν))U_ν(x-ν)
            minus *= dagger(
              neib.neighborLink(mu, nu, Direction::FORWARD, this->link_bc));
            minus *= dagger(
              it.neighborLink(nu, mu, Direction::BACKWARD, this->link_bc));
            minus *=
              it.neighborLink(nu, nu, Direction::BACKWARD, this->link_bc);

            A += plus.value() + minus.value();
          }
        }
        // MC part
        for (auto h = 0LU; h < hits; ++h) {
          auto U_new = update(U_old, gen, eps);
          auto delta_S = -beta * std::real((U_new.value() - U_old.value()) * A);
          auto boltzmann = std::exp(-1. * delta_S);

          auto change = this->accept(boltzmann, gen);

          // std::cout << delta_S << "\t" << boltzmann << "\t" << change <<
          // std::endl;

          if (change) {
            ++accepted;
            U_old = U_new;
            (*it)[mu] = U_new;
          }
        }
      }
    }

    // Return the number of accepted changes for acceptance tuning and
    // book keeping
    return accepted;
  }

  /// Overrelaxation
  void overrelaxation(size_t const &sweeps = 1) const
  {
    for (auto s = 0LU; s < sweeps; ++s) {
      // Iterate over lattice sites
      for (LinkLatticeIterator<LinkLattice<U1, dim>> it =
             this->lattice->begin();
           it != this->lattice->end(); ++it) {
        // Loop over links
        for (auto mu = 0LU; mu < this->lattice->dimensions(); ++mu) {
          // Compute staples
          std::complex<double> A = std::complex<double>(0.0);
          for (auto nu = 0LU; nu < this->lattice->dimensions(); ++nu) {
            if (nu != mu) {
              /// \todo Check boundary conditions in staples
              U1 plus(std::complex<double>(1., 0.));
              U1 minus(std::complex<double>(1., 0.));

              // U_ν(x+µ)dagger(U_µ(n+ν))dagger(U_ν(x))
              plus *=
                it.neighborLink(mu, nu, Direction::FORWARD, this->link_bc);
              plus *= dagger(
                it.neighborLink(nu, mu, Direction::FORWARD, this->link_bc));
              plus *= dagger(it.link(nu));

              // access to site (x-ν)
              auto neib = it.neighbor(nu, -1);

              // dagger(U_ν(x+µ-ν))dagger(U_µ(n-ν))U_ν(x-ν)
              minus *= dagger(
                neib.neighborLink(mu, nu, Direction::FORWARD, this->link_bc));
              minus *= dagger(
                it.neighborLink(nu, mu, Direction::BACKWARD, this->link_bc));
              minus *=
                it.neighborLink(nu, nu, Direction::BACKWARD, this->link_bc);

              A += plus.value() + minus.value();
            }
          }

          auto alpha = std::arg(A);
          auto phi = std::arg((it.link(mu)).value());

          // update phase
          phi = 2 * M_PI - 2 * alpha - phi;

          (*it)[mu] = std::polar(1., phi);
        }
      }
    }
  }

  template <typename generator>
  double Sweep(generator &gen, const double &eps = 0.5, size_t nsweep = 1,
               size_t nhit = 1, size_t nor = 1) const
  {

    auto accepted = 0LU;
    for (auto sw = 0ul; sw < nsweep; ++sw) {
      accepted += multihit_MC(gen, nhit, eps);
      overrelaxation(nor);
    }

    return static_cast<double>(accepted) /
           static_cast<double>(nhit * nsweep * lattice->volume() *
                               lattice->dimensions());
  }

  template <typename generator>
  void prepare_hot(generator &gen) const
  {

    // Iterate over lattice sites
    for (LinkLatticeIterator<LinkLattice<U1, dim>> it = this->lattice->begin();
         it != this->lattice->end(); ++it) {
      // Loop over links
      for (auto mu = 0LU; mu < this->lattice->dimensions(); ++mu) {
        (*it)[mu] = randomU1(gen);
      }
    }
  }
};


template <size_t dim>
class QEDGaugeFix {

 private:
  LinkLattice<U1, dim> *const lattice;
  double beta;


 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  /// Default constructor
  QEDGaugeFix() : beta(1.0) { lattice = nullptr; }

  /// Construct with given lattice
  QEDGaugeFix(LinkLattice<U1, dim> &orig, const double &coupling = 1.0)
    : lattice(&orig), beta(coupling)
  {
  }


  // **********************************************************************
  // Deconstructor
  // **********************************************************************

  ~QEDGaugeFix() = default;


  // **********************************************************************
  // Member Functions
  // **********************************************************************
  void randomGaugeTransformation()
  {
    auto &lat = *lattice;
    LinkLatticeIterator<LinkLattice<U1, dim>> lli(0, lat);

    // Init ranlux with true random number (Well, at least we try. The C++
    // standard does not guaranty std::random_device() gives a true random
    // number)
    std::random_device rd;
    std::seed_seq sseq({ rd(), rd(), rd() });
    // std::seed_seq sseq ({1,2,3,4});
    std::ranlux48 generator(sseq);


    // std::ostream_iterator<unsigned> out (std::cout," ");
    // sseq.param(out); std::cout << "\n" << generator << "\n" << std::endl;
    // Uniformly pick a phase from [-pi,pi)
    // Note that this is consistent with the std::arg function from std::complex
    std::uniform_real_distribution<double> dist(-1.0, 1.0);


    for (lli = lat.begin(); lli != lat.end(); ++lli) {
      double rphase = M_PI * dist(generator);

      std::complex<double> ii(0, 1);

      for (auto i = 0UL; i < dim; ++i) {
        auto n_dir = lli.neighbor(i, -1);
        (*lli)[i] *= std::exp(-ii * rphase);
        (*n_dir)[i] *= std::exp(ii * rphase);
      }
    }
  }


  double LandauGaugeFunctional() const
  {

    double res = 0;
    std::size_t i;
#pragma omp parallel for private(i) reduction(+ : res)
    for (auto idx = 0LU; idx < lattice->volume(); ++idx) {
      auto links = (*lattice)[idx];

      for (i = 0; i < lattice->dimensions(); ++i) {
        res += std::cos(links[i].phase());
      }
    }
    return res;
  }

  void localLandauGauge(size_t const &idx, const double or_param = 1.0)
  {
    auto &lat = *lattice;
    LinkLatticeIterator<LinkLattice<U1, dim>> lli(idx, lat);

    // Calculate the phase to locally fix the config to Landau gauge
    double num = 0.0;
    double den = 0.0;

    for (auto i = 0UL; i < dim; ++i) {
      auto n_dir = lli.neighbor(i, -1);
      auto here = (lli.link(i)).phase();
      auto neib = (n_dir.link(i)).phase();

      num += (std::sin(here) - std::sin(neib));
      den += (std::cos(here) + std::cos(neib));
    }

    auto phase = or_param * std::atan(num / den);

    // Apply gauge transformation to all adjacent links
    std::complex<double> ii(0, 1);
    for (auto i = 0UL; i < dim; ++i) {
      auto n_dir = lli.neighbor(i, -1);
      (*lli)[i] *= std::exp(-ii * phase);
      (*n_dir)[i] *= std::exp(ii * phase);
    }
  }

  double localLandauGaugeQuality(size_t const &idx)
  {
    auto &lat = *lattice;
    LinkLatticeIterator<LinkLattice<U1, dim>> lli(idx, lat);

    double res = 0.0;
    for (auto i = 0UL; i < dim; ++i) {
      auto n_dir = lli.neighbor(i, -1);
      auto here = (lli.link(i)).phase();
      auto neib = (n_dir.link(i)).phase();

      res += (std::sin(here) - std::sin(neib));
    }

    return res * res;
  }


  /// \Todo Implement boundary condition
  void LandauGaugeSweep(const double or_param = 1.0)
  {

    auto lat = *lattice;
    // FullLatticeIterator<FullLattice<CPtype,U1, dim>> fli(lat);

    // std::size_t count[] = {0,0};
    // Loop over all even and odd sites separately
    for (auto eo = 0; eo < 2; ++eo) {

#pragma omp parallel for
      for (auto idx = 0LU; idx < lat.volume(); ++idx) {
        auto coords = lat.linearIndexToCoord(idx);

        auto csum = std::accumulate(coords.begin(), coords.end(), 0);
        if (eo == csum % 2) {
          localLandauGauge(idx, or_param);
          // count[eo]++;
        }
      }
    }

    // std::cout << "Sites: " << count[0] << "/" << count[1] << std::endl;
    // std::cout << "Volume: " << lat.volume() << std::endl;
  }

  double LandauGaugeQuality()
  {

    auto lat = *lattice;

    double res = 0.0;
#pragma omp parallel for reduction(+ : res)
    for (auto idx = 0LU; idx < lat.volume(); ++idx) {
      res += localLandauGaugeQuality(idx);
    }

    return res / static_cast<double>(lat.volume());
  }

  std::size_t LandauGaugeDriver(const size_t gc_num, const size_t sw_num,
                                const double or_param = 1.0)
  {

    double minus_inf = std::numeric_limits<double>::lowest();
    LinkLattice<U1, dim> gauge_copy(*lattice), best_copy(*lattice);


    std::size_t iter = 0LU;

    QEDGaugeFix<dim> gc_gf(gauge_copy, this->beta);
    // double best= gc_gf.LandauGaugeQuality();
    // double last= std::numeric_limits<double>::max();

    double best = minus_inf; // Init with large negative number
    double last = 0;

    for (auto i = 0UL; i < gc_num; ++i) {


      if (0LU != i) // Use the original config once
      {
        gc_gf.randomGaugeTransformation();
      }


      // double last_LF=gc_gf.LandauGaugeFunctional();
      for (auto s = 0UL; s < sw_num; ++s) {
        gc_gf.LandauGaugeSweep(or_param);
        // last=gc_gf.LandauGaugeQuality();
        ++iter;

        if (gc_gf.LandauGaugeQuality() < 1.e-9) {
          // //last_LF = new_LF;
          // last=gc_gf.LandauGaugeQuality();
          last = gc_gf.LandauGaugeFunctional();
          // // std::cout << "\t\t Best local LF: " << last_LF << "\t (best
          // global: " << best << ")"
          // // 	  << "(" << ++s << " sweeps)" <<  std::endl;
          break;
        }
        // last_LF = new_LF;
      }

      // last = last_LF;
      if (last > best) {
        best_copy = gauge_copy;
        best = last;
      }
    }

    *lattice = best_copy;
    return iter;
  }
};
