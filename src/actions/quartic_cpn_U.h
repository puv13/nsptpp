/// \file
///
/// \brief Quartic CP(N) action with SU(N) fields as dynamical variables


#pragma once

#include "./cpn_common_aux.h"
#include "sun.h"

template <size_t N, size_t dim>
class QuarticCPNActionU {

 private:
  SiteLattice<nummat<N>, dim> const *lattice;
  double beta;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  /// Default constructor
  QuarticCPNActionU() : beta(1.0) { lattice = nullptr; }

  /// Construct with given lattice
  QuarticCPNActionU(SiteLattice<nummat<N>, dim> const &orig,
                    const double &coupling = 1.0)
    : lattice(&orig), beta(coupling)
  {
  }
}
