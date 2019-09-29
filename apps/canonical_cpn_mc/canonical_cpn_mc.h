///! \file
/// Header for canonical_cpn_mc.cpp

#pragma once

#include "../src/actions/canonical_cpn.h"
#include "../src/gaugegroups/u1.h"
#include "../src/lattice.h"
#include "../src/latticefields/cp.h"
#include "../src/progress_bar.h"
#include "../src/stat.h"
#include "./canonical_cpn_mc_params.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>


using size_t = std::size_t;
using clck = std::chrono::high_resolution_clock;


template <typename FullLatticeType, typename generator>
void canonical_CP_random_init(FullLatticeType &lattice, generator &gen);

template <typename T, size_t N, size_t dim, typename P = double>
class canonical_CP_update {

 private:
  FullLattice<CP<T, N>, U1, dim> *const lattice;
  double beta;
  BoundaryCondition<P, dim> site_bc;
  BoundaryCondition<double, dim> link_bc;


 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************
  canonical_CP_update()
    : lattice(nullptr), beta(1.0), site_bc(BoundaryCondition<P, dim>()),
      link_bc(BoundaryCondition<double, dim>())
  {
  }

  canonical_CP_update(
    FullLattice<CP<T, N>, U1, dim> *const l, double const &coupling,
    BoundaryCondition<P, dim> const &sb = BoundaryCondition<P, dim>(),
    BoundaryCondition<double, dim> const &lb = BoundaryCondition<double, dim>())
    : lattice(l), beta(coupling), site_bc(sb), link_bc(lb)
  {
  }


  // **********************************************************************
  // Destructor
  // **********************************************************************

  ~canonical_CP_update() = default;


  // **********************************************************************
  // Member functions
  // **********************************************************************

  /// Helper function for MC update
  template <typename generator>
  bool accept(double const &res, generator &gen)
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

  /// MC update of CP(N-1) field at given site
  template <typename generator>
  bool Site(size_t site, generator &gen, double const &eps = 1.e-2)
  {

    if (site > lattice->volume()) {
      throw std::runtime_error("Site out of range");
    }

    CanonicalCPNAction<std::complex<double>, N, dim> action(*lattice,
                                                            this->beta);
    // auto action_old = action.atSite(site);

    // Get iterator
    auto coord = lattice->linearIndexToCoord(site);
    auto it = lattice->at(coord);

    // Get CP(N-1) field  at site
    auto old_field = it.site();
    auto field = old_field;

    // MC Update
    field = update(field, gen, eps);

    // New action
    auto sum = std::real(T(0.0));
    auto diff = T(0.0);

    for (auto i = 0lu; i < dim; i++) {
      auto ilink = it.link(i).value();
      auto ilink_bwd = it.neighborLink(i, i, Direction::BACKWARD).value();

      auto field_fwd = it.neighborSite(i, Direction::FORWARD, this->site_bc);
      auto field_bwd = it.neighborSite(i, Direction::BACKWARD, this->site_bc);

      // Change in action due to change in n(x)^dagger
      diff = scalar_prod(field, field_fwd) - scalar_prod(old_field, field_fwd);
      sum += std::real(diff * ilink);

      // Change in action due to change in n(x)
      diff = scalar_prod(field_bwd, field) - scalar_prod(field_bwd, old_field);
      sum += std::real(diff * ilink_bwd);
    }

    auto action_diff = -2. * (this->beta) * N * sum;

    auto boltzmann = std::exp(-action_diff);

    bool change = this->accept(boltzmann, gen);


    if (change) {
      (*it).site = field;
    }
    return change;
  }

  /// MC update of U1 field at given site in given direction
  template <typename generator>
  bool Link(size_t site, size_t dir, generator &gen, double const &eps = 1.e-2)
  {

    if (site >= lattice->volume()) {
      throw std::runtime_error("Site out of range");
    }
    if (dir >= lattice->dimensions()) {
      throw std::runtime_error("Direction out of range");
    }

    // Get iterator
    auto coord = lattice->linearIndexToCoord(site);
    auto it = lattice->at(coord);

    // CanonicalCPNAction< std::complex<double>,N,dim>
    // action(*lattice,this->beta);
    // auto action_old = action.atSite(site);

    auto link_old = it.link(dir);
    auto link = update(link_old, gen, eps);

    auto n_field = it.neighborSite(dir, Direction::FORWARD, this->site_bc);
    auto field = it.site();

    // Change in action and weight
    auto diff_fast = std::real(scalar_prod(field, n_field) * link_old.value() -
                               scalar_prod(field, n_field) * link.value());
    auto boltzmann = std::exp(-2. * (this->beta) * N * diff_fast);

    bool change = this->accept(boltzmann, gen);

    if (change) {
      (*it).links[dir] = link;
    }

    return change;
  }

  template <typename generator>
  std::array<double, 2> Sweep(generator &gen, size_t nsweep = 1,
                              double const &seps = 5. / static_cast<double>(N),
                              double const &leps = 3.0)
  {

    size_t accepted = 0;
    size_t link_accepted = 0;
    for (auto sw = 0ul; sw < nsweep; ++sw) {

      for (auto site = 0ul; site < lattice->volume(); ++site) {

        // std::uniform_int_distribution<size_t> intdist(0,lattice->volume()-1);

        // Site update
        bool change = this->Site(site, gen, seps);
        if (change) {
          ++accepted;
        }

        // Link update
        for (auto nu = 0ul; nu < dim; ++nu) {
          change = this->Link(site, nu, gen, leps);
          if (change) {
            ++link_accepted;
          }
        }
      }
    }

    auto vol = lattice->volume();
    double s_updates = static_cast<double>(nsweep * vol);
    double l_updates = static_cast<double>(nsweep * vol * dim);

    double site_ratio = static_cast<double>(accepted) / s_updates;
    double link_ratio = static_cast<double>(link_accepted) / l_updates;

    return std::array<double, 2>{ site_ratio, link_ratio };
  }
};
