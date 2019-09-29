/// \file
///
/// \brief Principal chiral model

#pragma once

#include "../src/gaugegroups/sun.h"
#include "../src/lattice.h"
#include "action.h"
#include "expansion.h"
#include "nummat.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <numeric>

// Global constants known at compile time
// ----------------------------------------------------------------------
#ifdef SET_N_
constexpr std::size_t N = SET_N_;
#else
constexpr std::size_t N = 3;
#endif

// We need pi a lot here, so ...
#ifndef M_PI
#define M_PI 3.14159265358979323846 // Double precision pi
#endif

constexpr double PCM_Langevin_noise_scale = std::sqrt(2.);

using cmplx = std::complex<double>;

/// \brief Helper class for twisted PCM boundary conditions
template <size_t N>
class PCMtwistMatrix {
 private:
  nummat<N> bcmat;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  ///  \brief Standard Ctor
  ///  \details Set diagonal entries of bcmat to 1 <--> Periodic BC
  explicit PCMtwistMatrix() : bcmat(1.0) {}

  /// \brief Init Ctor
  /// \details Set diagonal entries of bcmat to given value
  PCMtwistMatrix(cmplx const &init) : bcmat(init) {}

  /// \brief Array Ctor
  /// \details Set bcmat according to a given matrix
  explicit PCMtwistMatrix(nummat<N> const &init_mat) : bcmat(init_mat) {}

  explicit PCMtwistMatrix(std::array<cmplx, N> const &init_ar)
    : bcmat(diag(init_ar))
  {
  }

  // **********************************************************************
  // Member Functions
  // **********************************************************************


  /// Access bcmat of PCMtwistMatrix
  nummat<N> get_mat() const
  {
    nummat<N> res(bcmat);
    return res;
  }

  void set_mat(nummat<N> &orig) { bcmat = orig; }

  /// Get size
  constexpr static size_t size() { return N; }


  /// Multiplication with PCM(N) field, i.e. with a nummat<N>, from the right
  nummat<N> operator*(nummat<N> const &field) const
  {

    nummat<N> res;

    // Multiply field from left with bcmat
    res = bcmat * field;

    // Multiply res from right with dagger(bcmat))
    res = res * bcmat.dagger();
    return res;
  }

  /// Multiplication with a scalar
  PCMtwistMatrix<N> operator*(cmplx const &rhs) const
  {
    auto res = bcmat;
    return rhs * res;
  }

  // **********************************************************************
  // Friend functions
  // **********************************************************************

  /// Multiplication with PCM(N) field, i.e. with a nummat<N>, from the
  /// left
  friend nummat<N> operator*(nummat<N> const &field,
                             PCMtwistMatrix<N> const &rhs)
  {

    nummat<N> res;

    // Multiply field from left with dagger(bcmat)
    res = rhs.bcmat.dagger() * field;

    // Multiply field from right with diag(phases)
    res = res * rhs.bcmat;

    return res;
  }
};


/// \brief Helper class for twisted PCM boundary conditions
template <typename T, size_t N>
class PCMtwistPhase {
 private:
  std::array<T, N> phases;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  ///  \brief Standard Ctor
  ///  \details Set all phases to 1 <--> Periodic BC
  explicit PCMtwistPhase() { phases.fill(T(1.0)); }

  /// \brief Init Ctor
  /// \details Set all phases to the given value
  explicit PCMtwistPhase(T const &init) { phases.fill(init); }

  /// \brief Array Ctor
  /// \details Set all phases according to a given array
  explicit PCMtwistPhase(std::array<T, N> const &init_ar)
  {
    for (size_t i = 0; i < N; ++i) {
      this->phases[i] = init_ar[i];
    }
  }

  // **********************************************************************
  // Member Functions
  // **********************************************************************


  /// Access phases of PCMtwistPhase
  T &operator[](std::size_t const idx) { return phases[idx]; }

  /// Access phases of constant PCMtwistPhase
  constexpr T const &operator[](const std::size_t idx) const
  {
    return phases[idx];
  }

  /// Get size
  constexpr static size_t size() { return N; }


  /// Multiplication with PCM(N) field, i.e. with a nummat<N>, from the right
  nummat<N> operator*(nummat<N> const &field) const
  {

    nummat<N> res;

    // Multiply field from left with diag(phases)
    for (auto i = 0LU; i < N; ++i) {
      for (auto j = 0LU; j < N; ++j) {
        res(i, j) = field(i, j) * this->phases[i];
      }
    }

    // Multiply field from right with dagger(diag(phases))
    for (auto i = 0LU; i < N; ++i) {
      for (auto j = 0LU; j < N; ++j) {
        res(i, j) *= std::conj(this->phases[j]);
      }
    }
    return res;
  }

  /// Multiplication with a scalar
  PCMtwistPhase<T, N> operator*(T const &rhs) const
  {
    auto res = *this;
    for (size_t i = 0; i < N; ++i) {
      res[i] *= rhs;
    }
    return res;
  }

  // **********************************************************************
  // Friend functions
  // **********************************************************************

  /// Multiplication with PCM(N) field, i.e. with a nummat<N>, from the
  /// left
  friend nummat<N> operator*(nummat<N> const &field,
                             PCMtwistPhase<T, N> const &rhs)
  {

    nummat<N> res;

    // Multiply field from left with dagger(diag(phases))
    for (auto i = 0LU; i < N; ++i) {
      for (auto j = 0LU; j < N; ++j) {
        res(i, j) = field(i, j) * std::conj(rhs[i]);
      }
    }

    // Multiply field from right with diag(phases)
    for (auto i = 0LU; i < N; ++i) {
      for (auto j = 0LU; j < N; ++j) {
        res(i, j) *= rhs[j];
      }
    }
    return res;
  }
};

template <typename T, size_t N, size_t dim>
class PCMAction {

 private:
  SiteLattice<nummat<N>, dim> const *lattice;
  double beta;
  typedef nummat<N> PCMtype;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  /// Default constructor
  PCMAction() : beta(1.0) { lattice = nullptr; }

  /// Construct with given lattice
  PCMAction(SiteLattice<nummat<N>, dim> const &orig,
            const double &coupling = 1.0)
    : lattice(&orig), beta(coupling)
  {
  }


  // **********************************************************************
  // Deconstructor
  // **********************************************************************
  ~PCMAction() = default;


  // **********************************************************************
  // Member Functions
  // **********************************************************************

  /// Setter for beta
  void setbeta(const double &coupling) { beta = coupling; }
  /// Getter for beta
  double getbeta() { return this->beta; }


  // Total Energy. The Definition is a bit wired, but follows PhysRevD.49.1621
  template <typename P = double>
  auto
  total(const BoundaryCondition<P, dim> &bc = BoundaryCondition<P, dim>()) const
    -> decltype(std::abs(T(1.)))
  {

    auto res = std::abs(T(0.0));

    // Loop over all sites
    for (auto fli = lattice->begin(); fli != lattice->end(); ++fli) {
      auto site = fli.site();

      for (auto i = 0lu; i < dim; i++) {
        auto sp = std::real(
          trace(site.mult_conj(fli.neighborSite(i, Direction::FORWARD, bc))));
        res += sp;
      }
    }

    return res;
  }


  // Energy density The Definition is a bit wired, but follows PhysRevD.49.1621
  template <typename P = double>
  auto energyDensity(
    const BoundaryCondition<P, dim> &bc = BoundaryCondition<P, dim>()) const
    -> decltype(std::abs(T(1.)))
  {
    auto vol = static_cast<double>(lattice->volume());
    return 1. -
           1. * total(bc) /
             (vol * static_cast<double>(dim) * static_cast<double>(N));
    // return 2.-1.*total(bc)/(vol*static_cast<double>(N));
  }
};


namespace PCMHelpers {

template <size_t N>
struct PCMCasimirStruct {
  const double site = static_cast<double>(N);
};


template <size_t N>
nummat<N> operator*(nummat<N> const &left, PCMCasimirStruct<N> const &CA)
{

  nummat<N> res;

  res = left * CA.site;

  return res;
}


template <size_t N>
nummat<N> operator*(PCMCasimirStruct<N> const &CA, nummat<N> const &right)
{
  nummat<N> res;
  res = CA.site * right;
  return res;
}
}


namespace PCMExpHelpers {


template <size_t N, size_t dim, typename Generator>
class randomPCMdf {
 public:
  nummat<N> operator()(Generator &rng) const
  {

    nummat<N> rd;
    rd = noiseSUN<N, Generator>(rng, PCM_Langevin_noise_scale);

    return rd;
  }
};


template <size_t N, size_t order>
auto operator*(Expansion<std::complex<double>, order> const &left,
               Expansion<nummat<N>, order> const &right)
  -> Expansion<nummat<N>, order>
{
  Expansion<nummat<N>, order> res;

  for (size_t n = 0; n < order; ++n) {
    for (size_t i = 0; i <= n; ++i) {
      res[n] += left[n - i] * right[i];
    }
  }
  return res;
}


template <size_t N, size_t order>
auto operator*(Expansion<nummat<N>, order> const &left,
               Expansion<std::complex<double>, order> const &right)
  -> Expansion<nummat<N>, order>
{
  return right * left;
}


template <size_t N, size_t order>
auto operator*(Expansion<std::complex<double>, order> const &left,
               nummat<N> const &right) -> Expansion<nummat<N>, order>
{
  Expansion<nummat<N>, order> res;

  for (size_t n = 0; n < order; ++n) {
    res[n] = left[n] * right;
  }

  return res;
}


template <size_t N, size_t order>
Expansion<nummat<N>, order> operator+(Expansion<nummat<N>, order> L,
                                      nummat<N> const &R)
{

  if (order > 1) {
    L[1] = L[1] + R;
  }

  return L;
}
}


// --------------------------------------------------------------------------------
// Drift Force
// --------------------------------------------------------------------------------
namespace NoExpansion {

using namespace PCMExpHelpers;

/// PCM Drift Force
template <size_t N, size_t dim, typename P = double>
class PCMForce
  : public Force<SiteLattice<nummat<N>, dim>, std::vector<nummat<N>>> {
 private:
  BoundaryCondition<P, dim> bc;

 public:
  typedef nummat<N> PCMdf;
  typedef std::vector<PCMdf> PCMForceField;
  typedef SiteLattice<nummat<N>, dim> PCMlattice;

  /// Returns force field  using bc boundary cond.
  PCMForceField operator()(PCMlattice const &lat) const
  {

    // Pre-allocate memory
    PCMForceField FF;
    FF.reserve(lat.volume());


    for (auto it = lat.begin(); it != lat.end(); ++it) {

      PCMdf res;

      const auto II = std::complex<double>(0.0, 1.0);

      auto unit_matrix = cmplx(0.0);
      for (auto i = 0lu; i < dim; i++) {
        auto field = it.site();

        auto field_fwd = it.neighborSite(i, Direction::FORWARD, bc);
        auto field_bwd = it.neighborSite(i, Direction::BACKWARD, bc);

        // First
        // ----------------------------------------
        auto aux = field.mult_conj(field_fwd);
        auto tr = trace(aux);
        res += aux;
        unit_matrix += tr;

        // Second (relative minus sign)
        // ----------------------------------------
        aux = field_bwd.mult_conj(field);
        tr = trace(aux);
        res -= aux;
        unit_matrix -= tr;
      }

      res *= (-1. * II * N);

      // std::cout << "Unit Matrix: " << unit_matrix << std::endl;
      res += ((unit_matrix)*nummat<N>(II));

      res += res.dagger();

      FF.push_back(0.5 * res);
      // ################################################################################
    }

    return FF;
  }

  /// Constructors
  explicit PCMForce() : bc(BoundaryCondition<P, dim>()) {}
  explicit PCMForce(const BoundaryCondition<P, dim> &_bc) : bc(_bc) {}
  ~PCMForce(){};
  /// Compiler complains about constexpr here ???
  bool isGroupValued() const { return false; }
};
}


namespace Expan {


/// PCM Drift Force
template <size_t N, size_t dim, size_t order, typename P = double>
class PCMForce : public Force<SiteLattice<Expansion<nummat<N>, order>, dim>,
                              std::vector<Expansion<nummat<N>, order>>> {
 private:
  BoundaryCondition<P, dim> bc;


 public:
  typedef Expansion<nummat<N>, order> PCMdf;
  typedef std::vector<PCMdf> PCMForceField;
  typedef SiteLattice<Expansion<nummat<N>, order>, dim> PCMlattice;


  /// Returns force field  using bc boundary cond.
  PCMForceField operator()(PCMlattice const &lat) const
  {

    // Pre-allocate memory
    PCMForceField FF;
    FF.reserve(lat.volume());


    for (auto it = lat.begin(); it != lat.end(); ++it) {

      PCMdf res;

      const auto II = std::complex<double>(0.0, 1.0);

      Expansion<std::complex<double>, order> unit_matrix(
        std::complex<double>(0.0));
      for (auto i = 0lu; i < dim; i++) {
        auto field = it.site();

        auto field_fwd = it.neighborSite(i, Direction::FORWARD, bc);
        auto field_bwd = it.neighborSite(i, Direction::BACKWARD, bc);

        // First
        // ----------------------------------------
        auto aux = field * dagger(field_fwd);
        auto tr = trace(aux);
        res += aux;
        unit_matrix += tr;

        // Second (relative minus sign)
        // ----------------------------------------
        aux = field_bwd * dagger(field);
        tr = trace(aux);
        res -= aux;
        unit_matrix -= tr;
      }

      res = res * (-1. * II * N);

      // std::cout << "Unit Matrix: " << unit_matrix << std::endl;
      // unit_matrix*=II;
      res += (nummat<N>(II) * unit_matrix);

      res += dagger(res);

      FF.push_back(0.5 * res);
      // ################################################################################
    }

    return FF;
  }

  /// Constructors
  explicit PCMForce() : bc(BoundaryCondition<P, dim>()) {}
  explicit PCMForce(const BoundaryCondition<P, dim> &_bc) : bc(_bc) {}
  ~PCMForce(){};
  /// Compiler complains about constexpr here ???
  bool isGroupValued() const { return false; }
};
}
