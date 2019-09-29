/// \file
///
/// Base classes for Action and Force
#pragma once
#include "expansion.h"
#include "lattice.h"

template <typename T, std::size_t order, std::size_t dimensions>
using Linkfield =
  std::array<SiteLattice<Expansion<T, order>, dimensions>, dimensions>;

/// Base class for actions
template <typename T, std::size_t order, std::size_t dimensions>
class Action {
 public:
  virtual ~Action(){};
  virtual double
  operator()(Linkfield<T, order, dimensions> const &gauge) const = 0;
};


/// Base class for drift forces
template <typename Latticetype, typename Forcefieldtype>
class Force {
 public:
  virtual ~Force(){};
  virtual Forcefieldtype operator()(Latticetype const &lat) const = 0;
  /// Should this be more general?
  constexpr bool isGroupValued() const;
};
